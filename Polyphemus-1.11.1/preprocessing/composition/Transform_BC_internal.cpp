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

int find_composition_id(Data<float,3>& total_mass, Data<float,4>& group_mass,
			  Data<float,3>& composition_bounds, int Nc ,int Ng, int z, int y, int x)
{
  //cout<<Nc<<"-"<<Ng<<">"<<t<<","<<z<<","<<y<<","<<x<<endl;
  float frac_group=0;
  Data<float,1> fgroup(Ng);
  for(int i = 0; i < Nc; i++)
  {
    int counter=0;
    if(Ng>1)
    {
      for(int g = 0; g < Ng; g++)
      {
	if(total_mass(z,y,x)>0.0)
	{
	  frac_group=group_mass(g, z, y,x)/total_mass(z,y,x);
	  fgroup(g)=frac_group;
	  //cout<<g<<" frac "<<frac_group<<" ["<<composition_bounds(i,g,0)
	  //<<","<<composition_bounds(i,g,1)<<"]"<<endl;
	}
	else
	  frac_group=0.0;
	if(composition_bounds(i,g,0)>0)
	{
	  if(frac_group>composition_bounds(i,g,0)&&frac_group<=composition_bounds(i,g,1))
	  {
	    counter++;
	  }
	}
	else
	{
	  if(frac_group>=composition_bounds(i,g,0)&&frac_group<=composition_bounds(i,g,1))
	    counter++;
	}
      }
      if(counter==Ng)
	return i;
    }
    else
      return 0;
  }
  cout<<"Error missing composition_bounds!"<<endl;
  cout<<" ("<<z<<","<<y<<","<<x<<")"<<endl;
  for(int g = 0; g < Ng; g++)
    cout<<"g("<<g<<")="<<fgroup(g)<<endl;
  abort();
  return Nc;
}


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
  int N_size_bins,Ncomposition,Nfrac,Ngroup;
  int Nt_out, Nz_out, Ny_out, Nx_out;
  real Delta_t_out, Delta_y_out, Delta_x_out,Delta_t_in;
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
  configuration.PeekValue("Delta_t_in", Delta_t_in);
  configuration.PeekValue("Delta_t_out", Delta_t_out);
  configuration.PeekValue("Delta_y", Delta_y_out);
  configuration.PeekValue("Delta_x", Delta_x_out);
  
  configuration.PeekValue("t_min", t_min_out);
  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);
  configuration.PeekValue("N_size_bins", N_size_bins);

  string vertical_levels_out;
  configuration.PeekValue("Vertical_levels", vertical_levels_out);

  // Input/output directories and files.
  string gas_file,aerosol_file; 
  string species_group_file;
  string composition_file,composition_data_file;
  string Directory_in, Directory_out;

  configuration.FindFromBeginning("[bc_files]");

  configuration.PeekValue("Species_Gas", gas_file);
  configuration.PeekValue("Species_Aerosol", aerosol_file);
  configuration.PeekValue("Species_Group", species_group_file);  
  configuration.PeekValue("Composition_configure", composition_file);//add composition_bounds configuration
  configuration.AddFile(composition_file);  
  configuration.PeekValue("Ncomposition", Ncomposition);
  configuration.PeekValue("N_frac", Nfrac);
  configuration.PeekValue("N_groups", Ngroup);  
  configuration.PeekValue("Composition_data", composition_data_file);
  Data<real,3> composition_bounds(Ncomposition,Ngroup,2);
  ExtStream composition_data(composition_data_file);
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
  RegularGrid<real> GridT_out(t_min_out, Delta_t_out, Nt_out);
  
  RegularGrid<real> GridZ_out(Nz_out);

  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  RegularGrid<real> GridComposition(Ncomposition);
  RegularGrid<real> GridGroup(Ngroup);
  // For boundary conditions.
  // All interfaces along z.
  RegularGrid<real> GridZ_all_interf_out(Nz_out + 1);
  // The interface where boundary concentration are required.
  RegularGrid<real> GridZ_interf_out(1);
  // Boundary layers along x and y.
  RegularGrid<real> GridY_interf_out(y_min_out - Delta_y_out / 2.,
				     (Ny_out+1) * Delta_y_out, 2);
  RegularGrid<real> GridX_interf_out(x_min_out - Delta_x_out / 2.,
				     (Nx_out+1) * Delta_x_out, 2);
  // Reads output altitudes.
  FormatText Heights_out;
  Heights_out.Read(vertical_levels_out, GridZ_all_interf_out);
  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_out(k) = (GridZ_all_interf_out(k) + GridZ_all_interf_out(k+1)) / 2.;
  // Sets the boundary layer altitude.
  GridZ_interf_out(0) = 1.5 * GridZ_all_interf_out(Nz_out-1)
    - 0.5 * GridZ_all_interf_out(Nz_out-2);
  // Sets time steps.??? ZS problem how to calculate
  
  //////////
  // DATA //
  //////////

  /////////////////////////
  // Boundary CONDITIONS //
  /////////////////////////

  int stating_t=int(t_min_out/Delta_t_in);
  for(int step=0;step<Nt_out;step++)
  {
    record=stating_t+int(step*Delta_t_out/Delta_t_in);
    cout<<"loop on time "<<record<<"s="<<step<<"<"<<Nt_out<<endl;
  // INITIAL GAS //
  ExtStream gas_species(gas_file);
  string gas_names;
// 
   FormatBinary<real> InputPolair;
// 
  if (!gas_species.is_open())
    throw string("Unable to open file \"") + gas_file + "\".";

  cout << "Gas:" << endl;

  while (gas_species.GetElement(gas_names))
    {
      // Concentrations and temporary arrays used to sum concentrations.
      Data<real, 3> Conc_out_x(GridZ_out, GridY_out, GridX_interf_out);
      Data<real, 3> Conc_out_y(GridZ_out, GridY_interf_out, GridX_out);
      Data<real, 3> Conc_out_z(GridZ_interf_out, GridY_out, GridX_out);
      Conc_out_x.SetZero();
      Conc_out_y.SetZero();
      Conc_out_z.SetZero();

      string File_conc = Directory_in + gas_names + "_x.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record,Conc_out_x);
      File_conc = Directory_in + gas_names + "_y.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record,Conc_out_y);
      File_conc = Directory_in + gas_names + "_z.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record,Conc_out_z);

      // Copy files without adding anything
      // Writes output files.
      FormatBinary<float> PolairOut;
      PolairOut.Append(Conc_out_x, Directory_out + gas_names+ "_x.bin");
      PolairOut.Append(Conc_out_y, Directory_out + gas_names+ "_y.bin");
      PolairOut.Append(Conc_out_z, Directory_out + gas_names+ "_z.bin");
    }

  gas_species.close();
  
  //INITIAL AEROSOL//
  
  for(int Nb=0;Nb<N_size_bins;Nb++)
  {
    ExtStream aerosol_species(aerosol_file);
    string aerosol_names;      
    if (!aerosol_species.is_open())
      throw string("Unable to open file \"") + aerosol_file + "\".";

    cout << "Aerosol size bin: " <<Nb<< endl;
    Data<real, 3> Conc_total_x(GridZ_out, GridY_out, GridX_interf_out);
    Data<real, 3> Conc_total_y(GridZ_out, GridY_interf_out, GridX_out);
    Data<real, 3> Conc_total_z(GridZ_interf_out, GridY_out, GridX_out);
    Conc_total_x.SetZero();
    Conc_total_y.SetZero();
    Conc_total_z.SetZero();
    Data<real, 4> Conc_group_x(GridGroup, GridZ_out, GridY_out, GridX_interf_out);
    Data<real, 4> Conc_group_y(GridGroup, GridZ_out, GridY_interf_out, GridX_out);
    Data<real, 4> Conc_group_z(GridGroup, GridZ_interf_out, GridY_out, GridX_out);
    Conc_group_x.SetZero();
    Conc_group_y.SetZero();
    Conc_group_z.SetZero();
    
    while (aerosol_species.GetElement(aerosol_names))
    {
      if(aerosol_names!="Number")
      {
      //cout << "   " << aerosol_names << endl;
      // input fields.
      // Concentrations and temporary arrays used to sum concentrations.
      Data<real, 3> Conc_out_x(GridZ_out, GridY_out, GridX_interf_out);
      Data<real, 3> Conc_out_y(GridZ_out, GridY_interf_out, GridX_out);
      Data<real, 3> Conc_out_z(GridZ_interf_out, GridY_out, GridX_out);
      
      char tostring[3];
      sprintf(tostring,"%d",Nb);
      string File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_x.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record,Conc_out_x);//problem here should use read record to speicified
      File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_y.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record, Conc_out_y);
      File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_z.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc,record, Conc_out_z);
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
      //compute the total bin/group mass for x
	  for (int z=0; z<Nz_out; z++)
	    for (int y=0; y<Ny_out; y++)
	      {	    
		Conc_total_x(z,y,0)+=Conc_out_x(z,y,0);
		Conc_group_x( group_id, z, y, 0)+=Conc_out_x(z,y,0);
		Conc_total_x(z,y,1)+=Conc_out_x(z,y,1);
		Conc_group_x(group_id, z, y, 1)+=Conc_out_x(z,y,1);
	      }
	  for (int z=0; z<Nz_out; z++)
	    for (int x=0; x<Nx_out; x++)
	      {
		Conc_total_y(z,0,x)+=Conc_out_y(z,0,x);
		Conc_group_y(group_id, z, 0, x)+=Conc_out_y(z,0,x);
		Conc_total_y(z,1,x)+=Conc_out_y(z,1,x);
		Conc_group_y(group_id, z, 1, x)+=Conc_out_y(z,1,x);
	      }
	  for (int y=0; y<Ny_out; y++)
	    for (int x=0; x<Nx_out; x++)
	      {
		Conc_total_z(0,y,x)+=Conc_out_z(0,y,x);
		Conc_group_z(group_id, 0, y, x)+=Conc_out_z(0,y,x);
	      }	      
      }
    }
    aerosol_species.close();
    
    //compute fraction of each group and locate composition_bounds of each grid
    real frac_group;
    Data<real, 1> fraction(Ngroup);
    Data<real, 3> Composition_grid_x(GridZ_out, GridY_out, GridX_interf_out);
    Data<real, 3> Composition_grid_y(GridZ_out, GridY_interf_out, GridX_out);
    Data<real, 3> Composition_grid_z(GridZ_interf_out, GridY_out, GridX_out);
    Composition_grid_x.SetZero();
    Composition_grid_y.SetZero();
    Composition_grid_z.SetZero();
	  for (int z=0; z<Nz_out; z++)
	    for (int y=0; y<Ny_out; y++)
	      {
		Composition_grid_x(z,y,0)=
		find_composition_id(Conc_total_x,Conc_group_x,
				    composition_bounds, Ncomposition , Ngroup,z,y,0);
		Composition_grid_x(z,y,1)=
		find_composition_id(Conc_total_x,Conc_group_x,
				    composition_bounds, Ncomposition , Ngroup,z,y,1);
	      }
	  for (int z=0; z<Nz_out; z++)
	    for (int x=0; x<Nx_out; x++)
	      {
		Composition_grid_y(z,0,x)=
		find_composition_id(Conc_total_y,Conc_group_y,
				    composition_bounds, Ncomposition , Ngroup,z,0,x);
		Composition_grid_y(z,1,x)=
		find_composition_id(Conc_total_y,Conc_group_y,
				    composition_bounds, Ncomposition , Ngroup,z,1,x);
	      }
	  for (int y=0; y<Ny_out; y++)
	    for (int x=0; x<Nx_out; x++)
	      {
		Composition_grid_z(0,y,x)=
		find_composition_id(Conc_total_z,Conc_group_z,
				    composition_bounds, Ncomposition , Ngroup,0,y,x);
	      }

    //redistribute 4d data into 5d based on composition_bounds id
    ExtStream aerosol_species2(aerosol_file);
    while (aerosol_species2.GetElement(aerosol_names))
    {
      // input fields.
      Data<real, 3> Conc_out_x(GridZ_out, GridY_out, GridX_interf_out);
      Data<real, 3> Conc_out_y(GridZ_out, GridY_interf_out, GridX_out);
      Data<real, 3> Conc_out_z(GridZ_interf_out, GridY_out, GridX_out);      
      // output fields.
      Data<real, 4> Conc_out_aer_x(GridComposition, GridZ_out, GridY_out, GridX_interf_out);
      Data<real, 4> Conc_out_aer_y(GridComposition, GridZ_out, GridY_interf_out, GridX_out);
      Data<real, 4> Conc_out_aer_z(GridComposition, GridZ_interf_out, GridY_out, GridX_out);
      Conc_out_aer_x.SetZero();
      Conc_out_aer_y.SetZero();
      Conc_out_aer_z.SetZero();
      char tostring[3];
      sprintf(tostring,"%d",Nb);
      string File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_x.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc, record,Conc_out_x);
      File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_y.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc, record, Conc_out_y);
      File_conc = Directory_in + aerosol_names+"_"+string(tostring)+ "_z.bin";
      //cout << "File_conc " << File_conc << endl;
      InputPolair.ReadRecord(File_conc, record, Conc_out_z);
      if(record==144&&aerosol_names=="Number")
	{
	cout<<aerosol_names<<endl;
	cout<<"x0=";
	Conc_out_x.PrintInfo();
	cout<<"y0=";
	Conc_out_y.PrintInfo();
	cout<<"z0=";
	Conc_out_z.PrintInfo();
	}
      int composition_id;
      //compute the total bin/group mass for x
	  for (int z=0; z<Nz_out; z++)
	    for (int y=0; y<Ny_out; y++)
	      {
		composition_id=Composition_grid_x(z,y,0);
		Conc_out_aer_x(composition_id,z,y,0)=Conc_out_x(z,y,0);
		composition_id=Composition_grid_x(z,y,1);
		Conc_out_aer_x(composition_id,z,y,1)=Conc_out_x(z,y,1);
	      }
	  for (int z=0; z<Nz_out; z++)
	    for (int x=0; x<Nx_out; x++)
	      {
		composition_id=Composition_grid_y(z,0,x);
		Conc_out_aer_y(composition_id,z,0,x)=Conc_out_y(z,0,x);
		composition_id=Composition_grid_y(z,1,x);
		Conc_out_aer_y(composition_id,z,1,x)=Conc_out_y(z,1,x);
	      }
	  for (int y=0; y<Ny_out; y++)
	    for (int x=0; x<Nx_out; x++)
	      {
		composition_id=Composition_grid_z(0,y,x);
		Conc_out_aer_z(composition_id,0, y, x)=Conc_out_z(0,y,x);
	      }
	// Output fields.
 	FormatBinary<float> PolairOut_aer;
	string File_out =Directory_out + aerosol_names +"_"+string(tostring)+ "_x.bin";
	PolairOut_aer.Append(Conc_out_aer_x,File_out );
	//cout << "File_out " << File_out << endl;
	File_out =Directory_out + aerosol_names +"_"+string(tostring)+ "_y.bin";
	PolairOut_aer.Append(Conc_out_aer_y,File_out );
	//cout << "File_out " << File_out << endl;
	File_out =Directory_out + aerosol_names +"_"+string(tostring)+ "_z.bin";
	PolairOut_aer.Append(Conc_out_aer_z,File_out );
	//cout << "File_out " << File_out << endl;
      if(record==144&&aerosol_names=="Number")
	{
	cout<<aerosol_names<<endl;
	cout<<"x1=";
	Conc_out_aer_x.PrintInfo();
	cout<<"y1=";
	Conc_out_aer_y.PrintInfo();
	cout<<"z1=";
	Conc_out_aer_z.PrintInfo();
	}
      }
    aerosol_species2.close();
   }
   record++;
}

  END;
  
  return 0;

}
