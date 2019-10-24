// Copyright (C) 2004-2005 CEREA
//     Authors: Marilyne Tombette and Denis Quélo
//
// CEREA (http://www.enpc.fr/cerea) is a joint laboratory of
// ENPC (http://www.enpc.fr) and EDF R&D (http://www.edf.fr).
//
// This file is part of a simulation system for air quality.
// 
// This code is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.


// Reads concentrations in Polair3D-like format on a coarse mesh (typically
// outputs from a previous simulation) and generates initial conditions for
// Polair3D on the finer mesh.


//////////////
// INCLUDES //

#include <iostream>

using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
using namespace AtmoData;
// INCLUDES //
//////////////


int main(int argc, char** argv)
{
 
  TRY;
  
  cout << endl;

  string main_config_file("Polair3D.cfg"), sec_config_file("");

  if (argc != 2 && argc != 3 && argc != 4 || !is_num(argv[argc - 1]))
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [main configuration file] [secondary config file] [date]\n";
      mesg += string("  ") + argv[0] + " [main configuration file] [date]\n";
      mesg += string("  ") + argv[0] + " [date]\n\n";
      mesg += "Arguments:\n";
      mesg += "  [main configuration file] (optional): main configuration file. Default: Polair3D.cfg\n";
      mesg += "  [secondary configuration file] (optional): secondary configuration file.\n";
      mesg += "  [date]: date in format YYYYMMDD.\n";
      cout << mesg << endl;
      return 1;
    }

  if (argc == 4)
    sec_config_file = argv[2];
  if (argc == 3 || argc == 4)
    main_config_file = argv[1];

  // Configuration files.
  if (!exists(main_config_file))
    throw string("Unable to find configuration file \"")
      + main_config_file + "\".";
  ConfigStreams configuration(main_config_file);
  if (exists(sec_config_file))
    configuration.AddFile(sec_config_file);


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////

  typedef float real;
  int date_mother_domain;
  // Output domain: Fine mesh (child domain).
  int N_size_section_aer,Ncomposition,Nfrac,Ngroup;
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_t_out, Delta_y_out, Delta_x_out;
  real t_min_out, y_min_out, x_min_out;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////

  cout << "Reading configuration files..."; cout.flush();

  configuration.FindFromBeginning("[domain]");
  configuration.PeekValue("Date", date_mother_domain);
  Date date_mother(date_mother_domain);
  Date date(convert<int>(argv[argc - 1]));  
  int record = date.GetDaysFrom(date_mother);
  configuration.PeekValue("Nt", Nt_out);
  configuration.PeekValue("Nz", Nz_out);
  configuration.PeekValue("Ny", Ny_out);
  configuration.PeekValue("Nx", Nx_out);

  configuration.PeekValue("Delta_t", Delta_t_out);
  configuration.PeekValue("Delta_y", Delta_y_out);
  configuration.PeekValue("Delta_x", Delta_x_out);
  
  configuration.PeekValue("t_min", t_min_out);
  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);
  configuration.PeekValue("N_size_section_aer", N_size_section_aer);

  string vertical_levels_out;
  configuration.PeekValue("Vertical_levels", vertical_levels_out);

  // Input/output directories and files.
  string gas_file,aerosol_file; 
  string species_group_file;
  string composition_file,composition_possibilities_file;
  string Directory_in, Directory_out;

  configuration.FindFromBeginning("[ic_files]");

  configuration.PeekValue("Species_Gas", gas_file);
  configuration.PeekValue("Species_Aerosol", aerosol_file);
  configuration.PeekValue("Species_Group", species_group_file);  
  configuration.PeekValue("Composition_configure", composition_file);//add composition configuration
  configuration.AddFile(composition_file);  
  configuration.PeekValue("Ncomposition", Ncomposition);
  configuration.PeekValue("N_frac", Nfrac);
  configuration.PeekValue("N_groups", Ngroup);  
  configuration.PeekValue("Composition_possibilities", composition_possibilities_file);
  Data<real,3> composition_bounds(Ncomposition,Ngroup,2);
  ExtStream composition_data(composition_possibilities_file);
  real number;
  int i=0;
  int j=Ngroup;
  int k=0;
  while(composition_data.GetNumber(number))
  {
    if(j==Ngroup)
    {
      i=number;
      j=0;
    }
    else
    {
      composition_bounds(i,j,k)=number;
      k++;
      if(k==2)
      {
	j++;
	k=0;
      }
    }
  }
  configuration.PeekValue("Database_Polair3D-ic", Directory_in);
  configuration.PeekValue("Directory_Polair3D-ic", Directory_out);
  cout << " Ncomposition=" <<Ncomposition<< endl;
  cout << " Nfrac="<< Nfrac<< endl;
  cout << " Ngroup=" <<Ngroup<< endl;  
  cout << " done." << endl;
  cout << endl;


  ///////////
  // GRIDS //
  ///////////
    
  cout << "Memory allocation for data fields..."; cout.flush();

  // Output grids.

  RegularGrid<real> GridZ_interf_out(Nz_out + 1);
  RegularGrid<real> GridZ_out(Nz_out);

  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  RegularGrid<real> GridComposition(Ncomposition);
  RegularGrid<real> GridGroup(Ngroup);
  ///////////////////////////
  // INPUT DATA PROCESSING //
  ///////////////////////////
  
  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels_out, GridZ_interf_out);
  for (int k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_interf_out(k) + GridZ_interf_out(k+1)) / 2.;

  /////////////////////////
  // INITIAL CONDITIONS //
  /////////////////////////
  // INITIAL GAS //
  ExtStream gas_species(gas_file);
  string gas_names;

  FormatBinary<real> InputPolair;

  if (!gas_species.is_open())
    throw string("Unable to open file \"") + gas_file + "\".";

  cout << "Gas:" << endl;

  while (gas_species.GetElement(gas_names))
    {
      // input fields.
      // Initial Conditions
      Data<real, 3> Conc_ic_out(GridZ_out, GridY_out,GridX_out);
      
      string File_conc = Directory_in + gas_names + ".bin";
      cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc, record,Conc_ic_out);

	  // Writes output files.
	  FormatBinary<float> PolairOut;
	  PolairOut.Append(Conc_ic_out, Directory_out + gas_names 
		       + ".bin");
    }

  gas_species.close();
  
  //INITIAL AEROSOL//
  
  for(int Nb=0;Nb<N_size_section_aer;Nb++)
  {
    ExtStream aerosol_species(aerosol_file);
    string aerosol_names;      
    if (!aerosol_species.is_open())
      throw string("Unable to open file \"") + aerosol_file + "\".";

    cout << "Aerosol:" << endl;
    Data<real, 3> Conc_total(GridZ_out, GridY_out,GridX_out);
    Conc_total.SetZero();
    Data<real, 4> Conc_group(GridGroup,GridZ_out, GridY_out,GridX_out);
    Conc_group.SetZero();
   
    while (aerosol_species.GetElement(aerosol_names))
    {
      //cout << "   " << aerosol_names << endl;
      // input fields.
      // Initial Conditions
      Data<real, 3> Conc_ic_out(GridZ_out, GridY_out,GridX_out);
      Data<real, 4> Conc_ic_out_compos(GridComposition,GridZ_out, GridY_out
	  ,GridX_out);
      Conc_ic_out_compos.SetZero();      
      char tostring[3];
      sprintf(tostring,"%d",Nb);
      string File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ ".bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc, record, Conc_ic_out);
      if(aerosol_names!="Number")
      {
      //find group_id based on species
      ExtStream groups(species_group_file);
      string groups_names;
      int group_id;
      while (groups.GetElement(groups_names))
      {
	if(aerosol_names==split(groups_names, "-")[0])
	{
	  group_id=convert<int>(split(groups_names, "_")[1]);
	  group_id=group_id-1;
	}
      }
      //compute the total bin/group mass
      for(int z = 0; z < Nz_out; z++)
	 for(int y = 0; y < Ny_out; y++)
	    for(int x = 0; x < Nx_out; x++)
	    {
	      if(Conc_ic_out(z,y,x)>0.0)
	      {
	      Conc_total(z,y,x)+=Conc_ic_out(z,y,x);
	      Conc_group(group_id,z, y,x)+=Conc_ic_out(z,y,x);
// 	      cout<<Conc_group(group_id,z, y,x)<<"/"<<Conc_total(z,y,x)<<"+"<<Conc_ic_out(z,y,x)
// 	         <<" ("<<z<<","<<y<<","<<x<<")"<<endl;
	      }
	    }
	//redistribute 3d data into 4d based on composition id
	// Initialise Conc_ic_out_compos from Conc_ic_out
	int Composition_id;
	if(group_id==Ngroup-1)
	{
	  Composition_id=0;
	}
	else
	{
	  for(int i=1;i<Ncomposition;i++)
	  {
	    if(composition_bounds(i,group_id,1)==1)//looking for pure external composition
	      Composition_id=i;
	  }
	}
	//cout<<"Composition_id "<<Composition_id<<endl;
	for(int i=0;i<Ncomposition;i++)
	{
	  if(i==Composition_id)
	  {
	    for(int z = 0; z < Nz_out; z++)
	      for(int y = 0; y < Ny_out; y++)
		for(int x = 0; x < Nx_out; x++)
		{
		  //cout<<GridComposition(i)<<" "<<GridZ_out(z)<<" "<<GridY_out(y)
		  //<<" "<<GridX_out(x)<<endl;
		  Conc_ic_out_compos(i,z,y,x)=Conc_ic_out(z,y,x);
		}
	    }
	  }
    }
    else//redistribute number concentrations based on group mass fractions
    {
	// Initialise Conc_ic_out_compos from Conc_ic_out      
      for(int z = 0; z < Nz_out; z++)
	for(int y = 0; y < Ny_out; y++)
	  for(int x = 0; x < Nx_out; x++)
	  {
	    for(int g = 0; g < Ngroup; g++)
	    {
	      real frac_group=0.0;
	      if(Conc_total(z,y,x)>0.0)
	      {
		frac_group=Conc_group(g,z, y,x)/Conc_total(z,y,x);
	      }
	      else
		frac_group=0.0;
	      int Composition_id;
	      if(g==Ngroup-1)
	      {
		Composition_id=0;
	      }
	      else
	      {
		for(int i=1;i<Ncomposition;i++)
		{
		  if(composition_bounds(i,g,1)==1)//looking for pure external composition
		  Composition_id=i;
		}
	      }
	      Conc_ic_out_compos(Composition_id,z,y,x)=frac_group*Conc_ic_out(z,y,x);
	    }
	  }
    }
	  // Writes output files.
 	  FormatBinary<float> PolairOut_aer;
	  string File_out =Directory_out + aerosol_names +"_"+string(tostring)+ ".bin";
 	  PolairOut_aer.Append(Conc_ic_out_compos,File_out );
	  cout << "File_out " << File_out << endl;	    
   }
  aerosol_species.close();
  }
  cout << endl;  
  END;
  
  return 0;

}
