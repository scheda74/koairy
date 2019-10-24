// Copyright (C) 2006-2008, ENPC - INRIA - EDF R&D
// Author(s): Marilyne Tombette
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// Polyphemus is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// Polyphemus is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/


//////////////
// INCLUDES //


#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4

#include "SeldonData.hxx"
using namespace SeldonData;

#include "AtmoData.hxx"
using namespace AtmoData;

#include "Optics.cxx"


// INCLUDES //
//////////////


template<class T>
void abs_(T& x)
{
  if (x < T(0))
    x = -x;
}


int main(int argc, char** argv)
{

  TRY;

  cout << endl;

  string main_config_file("optics.cfg"), sec_config_file("");

  if (argc != 4 && argc != 3)
    {
      string mesg = "Usage:\n";
      mesg += string("  ") + string(argv[0]);
      mesg += " [main configuration file] [secondary config file] [date]\n";
      mesg += string("  ") + string(argv[0]);
      mesg += " [main configuration file] [date]\n";
      mesg += string("  ") + string(argv[0]) + " [date]\n\n";
      mesg += "Arguments:\n";
      mesg += "  [main configuration file] (optional): ";
      mesg += "main configuration file. Default: optics.cfg\n";
      mesg += "  [secondary configuration file] (optional): ";
      mesg += "secondary configuration file.\n";
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

  Date date(convert<int>(argv[argc - 1]));

  cout << "Date: " << date.GetDate() << endl << endl;


  ////////////////////////
  // FIRST DECLARATIONS //
  ////////////////////////


  typedef float real;

  const real dry_aerosol_density = 1.4e6; // g.m^{-3}.

  int h, i, j, k, l, t, is, is2;

  // Input and output directories and files.
  string directory_simulation, directory_result;
  string file_temperature, file_pressure, file_specific_humidity;
  string directory_OPAC, file_water_refractive_index, file_species;
  string directory_efficiency_factor;

  // Input dimensions for meteological data.
  int Nt_in, Nz_in, Ny_in, Nx_in;

  // Input-data steps.
  real Delta_t_in, Delta_y_in, Delta_x_in;
  real t_min_in, y_min_in, x_min_in;

  // Simulation dimensions.
  int Nt_out, Nz_out, Ny_out, Nx_out;

  // Output steps.
  real Delta_t_out, Delta_y_out, Delta_x_out;
  real y_min_out, x_min_out;

  // Output altitudes.
  string file_level;

  // Refractive index of black carbon.
  int black_carbon_index(-999);

  // If the number of bins used in simulation is less than the number of bins
  // used for interpolation, then concentrations are interpolated on
  // Nbins_for_interpolation bins.
  int Nbins_for_interpolation = 10;

  int Nbins_in_simulation;


  /////////////////////////
  // CONFIGURATION FILES //
  /////////////////////////


  cout << "Reading configuration files...";
  cout.flush();

  configuration.Find("[domain]");

  // Dates.
  string date_meteo_str;
  configuration.PeekValue("Date", date_meteo_str);
  Date date_meteo(date_meteo_str);

  int record = date.GetDaysFrom(date_meteo);

  configuration.PeekValue("Nt", Nt_in);
  configuration.PeekValue("Nz", Nz_in);
  configuration.PeekValue("Ny", Ny_in);
  configuration.PeekValue("Nx", Nx_in);

  configuration.PeekValue("Delta_t", Delta_t_in);
  configuration.PeekValue("Delta_y", Delta_y_in);
  configuration.PeekValue("Delta_x", Delta_x_in);

  configuration.PeekValue("t_min", t_min_in);
  configuration.PeekValue("y_min", y_min_in);
  configuration.PeekValue("x_min", x_min_in);
  configuration.PeekValue("Vertical_levels", file_level);

  // Input/output files.
  configuration.FindFromBeginning("[paths]");

  configuration.PeekValue("Directory_simulation_result",
                          directory_simulation);

  configuration.PeekValue("File_temperature", file_temperature);
  configuration.PeekValue("File_pressure", file_pressure);
  configuration.PeekValue("File_specific_humidity", file_specific_humidity);

  configuration.PeekValue("Directory_efficiency_factor",
                          directory_efficiency_factor);

  configuration.PeekValue("Directory_OPAC", directory_OPAC);
  configuration.PeekValue("File_index_water", file_water_refractive_index);
  configuration.PeekValue("File_species_match", file_species);

  configuration.PeekValue("Directory_result", directory_result);

  configuration.Find("[domain_result]");

  // Dates.
  string date_result_str;
  configuration.PeekValue("Date", date_result_str);
  Date date_result(date_result_str);

  int record_result = date.GetDaysFrom(date_result);

  configuration.PeekValue("Nt", Nt_out);
  configuration.PeekValue("Nz", Nz_out);
  configuration.PeekValue("Ny", Ny_out);
  configuration.PeekValue("Nx", Nx_out);

  configuration.PeekValue("Delta_t", Delta_t_out);
  configuration.PeekValue("Delta_y", Delta_y_out);
  configuration.PeekValue("Delta_x", Delta_x_out);

  configuration.PeekValue("y_min", y_min_out);
  configuration.PeekValue("x_min", x_min_out);

  // Optics.
  vector<real> wavelength;
  int Nwavelength;
  int tabulation_index_real;
  int tabulation_index_imaginary;
  int Ndiameter, N_OPAC_wavelength, Nwater_wavelength;

  configuration.FindFromBeginning("[optic]");

  // Wavelengths of aerosol optical thickness.
  configuration.Find("Wavelength");
  split(configuration.GetLine(), wavelength);
  Nwavelength = wavelength.size();

  configuration.PeekValue("Tabulation_refractive_index_real",
                          tabulation_index_real);
  configuration.PeekValue("Tabulation_refractive_index_imaginary",
                          tabulation_index_imaginary);
  configuration.PeekValue("Ndiameter", Ndiameter);

  configuration.PeekValue("N_OPAC_wavelength", N_OPAC_wavelength);
  configuration.PeekValue("N_water_wavelength", Nwater_wavelength);

  // Aerosol.
  int Nbins, Nspecies;
  real min_diameter, max_diameter;
  string aerosol_water_name;

  configuration.FindFromBeginning("[aerosol]");
  configuration.PeekValue("Nbins", Nbins);
  configuration.PeekValue("Min_diameter", min_diameter);
  configuration.PeekValue("Max_diameter", max_diameter);
  configuration.PeekValue("Aerosol_water_name", aerosol_water_name);
  configuration.PeekValue("Nspecies", Nspecies);

  if (Nbins >= Nbins_for_interpolation)
    Nbins_in_simulation = Nbins;
  else
    Nbins_in_simulation = Nbins_for_interpolation;

  // Options.
  int option_dry_diameter, option_wet_index;
  int option_well_mixed_index, option_black_carbon_treatment;

  configuration.Find("[option]");

  configuration.PeekValue("Dry_diameter_option", option_dry_diameter);
  configuration.PeekValue("Wet_computation_option", option_wet_index);
  configuration.PeekValue("Well_mixed_computation_option",
                          option_well_mixed_index);
  configuration.PeekValue("Black_carbon_treatment",
                          option_black_carbon_treatment);

  cout << " done." << endl;
  cout << endl;


  ///////////
  // GRIDS //
  ///////////


  cout << "Memory allocation for data fields...";
  cout.flush();

  /*** Input settings. ***/

  // Input grids.
  RegularGrid<real> GridT_in(t_min_in, Delta_t_in, Nt_in);
  RegularGrid<real> GridZ_in(Nz_in);
  RegularGrid<real> GridY_in(y_min_in, Delta_y_in, Ny_in);
  RegularGrid<real> GridX_in(x_min_in, Delta_x_in, Nx_in);

  /*** Output settings. ***/

  // Output grids.
  RegularGrid<real> GridT_out(0.0, Delta_t_out, Nt_out);
  RegularGrid<real> GridZ_out(Nz_out);
  RegularGrid<real> GridZ_interf_out(Nz_out + 1);
  RegularGrid<real> GridY_out(y_min_out, Delta_y_out, Ny_out);
  RegularGrid<real> GridX_out(x_min_out, Delta_x_out, Nx_out);

  RegularGrid<real> GridSpecies(Nspecies);
  //  Without black carbon.
  RegularGrid<real> GridSpeciesNoBlackCarbon(Nspecies - 1);

  RegularGrid<real> GridBins(Nbins);
  RegularGrid<real> GridBins_in_simulation(Nbins_in_simulation);

  // Reads output altitudes.
  FormatText height_out;
  height_out.Read(file_level.c_str(), GridZ_interf_out);

  // Sets values at nodes.
  for (k = 0; k < Nz_out; k++)
    GridZ_in(k) = GridZ_out(k) = (GridZ_interf_out(k)
                                  + GridZ_interf_out(k + 1)) / 2.0;

  RegularGrid<real> GridWavelength(Nwavelength);
  for (i = 0; i < Nwavelength; i++)
    GridWavelength(i) = wavelength[i];


  //////////
  // DATA //
  //////////


  // Input fields.
  Data<real, 4> TemperatureIn(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> PressureIn(GridT_in, GridZ_in, GridY_in, GridX_in);
  Data<real, 4> SpecificHumidityIn(GridT_in, GridZ_in, GridY_in, GridX_in);

  Data<real, 4> Temperature(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> Pressure(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> SpecificHumidity(GridT_out, GridZ_out, GridY_out, GridX_out);
  Data<real, 4> RelativeHumidity(GridT_out, GridZ_out, GridY_out, GridX_out);

  cout << " done." << endl;
  cout << endl;


  /////////////////
  // READS INPUT //
  /////////////////


  cout << "Extracting meteo data...";
  cout.flush();

  FormatBinary<float> Input;

  Input.ReadRecord(file_temperature, record, TemperatureIn);
  Input.ReadRecord(file_pressure, record, PressureIn);
  Input.ReadRecord(file_specific_humidity, record, SpecificHumidityIn);

  cout << " done." << endl;

  cout << "Interpolating meteo data on simulation grid...";
  cout.flush();
  LinearInterpolationRegular(TemperatureIn, Temperature);
  LinearInterpolationRegular(PressureIn, Pressure);
  LinearInterpolationRegular(SpecificHumidityIn, SpecificHumidity);
  cout << " done." << endl;

  cout << "Computing relative humidity...";
  cout.flush();
  ComputeRelativeHumidity(SpecificHumidity, Temperature, Pressure,
                          RelativeHumidity);
  cout << " done." << endl;

  cout << endl;

  cout << "Reading input data files...";
  cout.flush();

  // Refractive index for pure species when relative humidity is zero.
  Data<string, 1, real> SpeciesNames(GridSpecies);
  Data<string, 1, real> OPACNames(GridSpecies);
  RegularGrid<real> GridWavelengthOPAC(N_OPAC_wavelength);
  Data<real, 2> IndexReal(GridSpecies, GridWavelengthOPAC);
  Data<real, 2> IndexImaginary(GridSpecies, GridWavelengthOPAC);

  FormatFormattedText match_OPAC(string("<e><e>"));
  match_OPAC.SetDelimiters("\t");
  match_OPAC.Read(file_species, "0", SpeciesNames);
  match_OPAC.Read(file_species, "1", OPACNames);
  for (i = 0; i < Nspecies; i++)
    {
      if (SpeciesNames(i) == "PBC")
        black_carbon_index = i;
      ExtStream  file_OPAC(directory_OPAC + string("/") + OPACNames(i));
      for (k = 0; k < 17; k++)
        file_OPAC.GetFullLine();
      for (j = 0; j < N_OPAC_wavelength; j++)
        {
          file_OPAC.GetElement(GridWavelengthOPAC(j));
          file_OPAC.SkipElements(6);
          file_OPAC.GetElement(IndexReal(i, j));
          file_OPAC.GetElement(IndexImaginary(i, j));
        }
      file_OPAC.Close();
    }

  // Refractive index for pure species when relative humidity is zero.
  string description;
  RegularGrid<real> GridWaterWavelength(Nwater_wavelength);
  Data<real, 1> WaterIndexInputReal(GridWaterWavelength);
  Data<real, 1> WaterIndexInputImaginary(GridWaterWavelength);
  FormatText input_tab(",");
  ifstream water_refractive_index_data(file_water_refractive_index.c_str());
  getline(water_refractive_index_data, description);
  input_tab.Read(water_refractive_index_data, GridWaterWavelength);
  input_tab.Read(water_refractive_index_data, WaterIndexInputReal);
  input_tab.Read(water_refractive_index_data, WaterIndexInputImaginary);
  water_refractive_index_data.close();

  // Tabulation of efficiency factors.
  string file_efficiency_factor, file_grid_Mie;
  RegularGrid<real> GridIndexReal(tabulation_index_real);
  RegularGrid<real> GridIndexImaginary(tabulation_index_imaginary);
  RegularGrid<real> GridDiameter(Ndiameter);

  file_grid_Mie = directory_efficiency_factor + "/Grid_Mie.dat";
  ifstream stream_efficiency_grid(file_grid_Mie.c_str());
  stream_efficiency_grid >> description;
  input_tab.Read(stream_efficiency_grid, GridIndexReal);
  input_tab.Read(stream_efficiency_grid, GridIndexImaginary);
  input_tab.Read(stream_efficiency_grid, GridDiameter);
  stream_efficiency_grid.close();
  cout << " done." << endl;

  cout << "Computing dry diameters...";
  cout.flush();

  Data<real, 1> dry_diameter(Nbins), dry_diameter_bound(Nbins + 1);
  real delta_diameter = log(max_diameter / min_diameter) / real(Nbins);

  for (i = 0; i < Nbins + 1; i++)
    dry_diameter_bound(i) = min_diameter * exp(real(i) * delta_diameter);

  for (i = 0; i < Nbins; i++)
    dry_diameter(i) = sqrt(dry_diameter_bound(i) * dry_diameter_bound(i + 1));

  Data<real, 1> dry_diameter_computed(Nbins_in_simulation);
  Data<real, 1> dry_diameter_bound_computed(Nbins_in_simulation + 1);
  real delta_diameter_computed = log(max_diameter / min_diameter)
    / real(Nbins_in_simulation);

  for (i = 0; i < Nbins_in_simulation + 1; i++)
    dry_diameter_bound_computed(i) = min_diameter
      * exp(real(i) * delta_diameter_computed);

  for (i = 0; i < Nbins_in_simulation; i++)
    dry_diameter_computed(i) = sqrt(dry_diameter_bound_computed(i)
                                    * dry_diameter_bound_computed(i + 1));

  Data<real, 2> coefficient_computed(Nbins_in_simulation, Nbins);
  coefficient_computed.Fill(0.0);
  for (i = 0; i < Nbins_in_simulation; i++)
    for (is = 0; is < Nbins; is++)
      if (dry_diameter_bound(is) <= dry_diameter_bound_computed(i))
        {
          if (dry_diameter_bound(is + 1) > dry_diameter_bound_computed(i))
            if (dry_diameter_bound(is + 1)
                <= dry_diameter_bound_computed(i + 1))
              coefficient_computed(i, is)
                += (dry_diameter_bound(is + 1)
                    - dry_diameter_bound_computed(i))
                / (dry_diameter_bound(is + 1) - dry_diameter_bound(is));
            else
              coefficient_computed(i, is)
                += (dry_diameter_bound_computed(i + 1)
                    - dry_diameter_bound_computed(i))
                / (dry_diameter_bound(is + 1) - dry_diameter_bound(is));
        }
      else if (dry_diameter_bound(is) < dry_diameter_bound_computed(i + 1))
        if (dry_diameter_bound(is + 1)
            < dry_diameter_bound_computed(i + 1))
          coefficient_computed(i, is) = 1.;
        else
          coefficient_computed(i, is) +=
            (dry_diameter_bound_computed(i + 1)
             - dry_diameter_bound(is))
            / (dry_diameter_bound(is + 1) - dry_diameter_bound(is));
  cout << " done." << endl;


  /////////////
  // INDICES //
  /////////////


  cout << "Reading aerosol concentrations...";
  cout.flush();

  Data<real, 6> ConcentrationAerosol(GridSpecies, GridBins_in_simulation,
                                     GridT_out, GridZ_out,
                                     GridY_out, GridX_out);
  ConcentrationAerosol.Fill(0.0);
  Data<real, 5> ConcentrationWater(GridBins_in_simulation,
                                   GridT_out, GridZ_out,
                                   GridY_out, GridX_out);
  ConcentrationWater.Fill(0.0);
  string file_concentration;

  Data<real, 4> concentration_tmp(GridT_out, GridZ_out, GridY_out, GridX_out);

  for (l = 0; l < Nspecies; l++)
    for (is = 0; is < Nbins; is++)
      {
        file_concentration = directory_simulation + string("/")
          + SpeciesNames(l) + string("_") + to_str(is) + ".bin";

        Input.ReadRecord(file_concentration, record_result,
                         concentration_tmp);
        for (t = 0; t < Nt_out; t++)
          for (h = 0; h < Nz_out; h++)
            for (j = 0; j < Ny_out; j++)
              for (i = 0; i < Nx_out; i++)
                for (is2 = 0; is2 < Nbins_in_simulation; is2++)
                  ConcentrationAerosol(l, is2, t, h, j, i)
                    += coefficient_computed(is2, is)
                    * concentration_tmp(t, h, j, i);
      }

  if (option_dry_diameter == 3 || option_wet_index == 2)
    for (is = 0; is < Nbins; is++)
      {
        file_concentration = directory_simulation + string("/")
          + aerosol_water_name + string("_") + to_str(is) + ".bin";
        Input.ReadRecord(file_concentration, record_result,
                         concentration_tmp);
        for (t = 0; t < Nt_out; t++)
          for (h = 0; h < Nz_out; h++)
            for (j = 0; j < Ny_out; j++)
              for (i = 0; i < Nx_out; i++)
                for (is2 = 0; is2 < Nbins_in_simulation; is2++)
                  ConcentrationWater(is2, t, h, j, i) +=
                    coefficient_computed(is2, is)
                    * concentration_tmp(t, h, j, i);
      }

  cout << "done." << endl;

  // Interpolates refractive index of pure species on desired wavelengths.
  Data<real, 2> PureSpeciesIndexReal(GridSpecies, GridWavelength);
  Data<real, 2> PureSpeciesIndexImaginary(GridSpecies, GridWavelength);

  for (k = 0; k < Nwavelength; k++)
    {
      int index_in_OPAC_table(-999);
      real coeff1(-999.), coeff2(-999.);
      for (i = 0; i < N_OPAC_wavelength - 1; i++)
        if (GridWavelengthOPAC(i) <= GridWavelength(k)
            && GridWavelengthOPAC(i + 1) > GridWavelength(k))
          {
            index_in_OPAC_table = i;
            coeff2 = (GridWavelength(k) - GridWavelengthOPAC(i))
              / (GridWavelengthOPAC(i + 1) - GridWavelengthOPAC(i));
            coeff1 = (GridWavelengthOPAC(i + 1) - GridWavelength(k))
              / (GridWavelengthOPAC(i + 1) - GridWavelengthOPAC(i));

            for (j = 0; j < Nspecies; j++)
              {
                PureSpeciesIndexReal(j, k) =
                  coeff1 * IndexReal(j, index_in_OPAC_table)
                  + coeff2 * IndexReal(j, index_in_OPAC_table + 1);

                PureSpeciesIndexImaginary(j, k) =
                  coeff1 * IndexImaginary(j, index_in_OPAC_table)
                  + coeff2 * IndexImaginary(j, index_in_OPAC_table + 1);
              }
          }

      if (index_in_OPAC_table == -999)
        throw string("Unable to find required wavelength in OPAC.");
    }

  Data<real, 1> WaterIndexReal(GridWavelength);
  Data<real, 1> WaterIndexImaginary(GridWavelength);
  LinearInterpolationRegular(WaterIndexInputReal, WaterIndexReal);
  LinearInterpolationRegular(WaterIndexInputImaginary, WaterIndexImaginary);

  // For Mie calculations, imaginary indices are positive.
  PureSpeciesIndexImaginary.Apply(abs_);

  real wet_diameter, wet_diameter_without_black_carbon;

  Data<real, 4> OpticalThickness(GridWavelength, GridT_out, GridY_out,
                                 GridX_out);
  Data<real, 4> OpticalThicknessScattering(GridWavelength, GridT_out,
                                           GridY_out, GridX_out);
  Data<real, 4> SingleScatteringAlbedo(GridWavelength, GridT_out, GridY_out,
                                       GridX_out);
  Data<real, 5> ExtinctionCoefficient(GridWavelength, GridT_out, GridZ_out,
                                      GridY_out, GridX_out);
  Data<real, 5> AbsorptionCoefficient(GridWavelength, GridT_out, GridZ_out,
                                      GridY_out, GridX_out);
  OpticalThickness.SetZero();
  for (k = 0; k < Nwavelength; k++)
    {
      cout << "Computing optical thickness at " << GridWavelength(k) * 1.e3
           << " nm..." << endl;
      cout << "    Reading efficiency factors file...";
      cout.flush();
      Data<real, 3> AbsorptionEfficiencyFactor(GridIndexReal,
                                               GridIndexImaginary,
                                               GridDiameter);
      Data<real, 3> ExtinctionEfficiencyFactor(GridIndexReal,
                                               GridIndexImaginary,
                                               GridDiameter);

      file_efficiency_factor = directory_efficiency_factor
        + string("/efficiency_factors_tab_")
        + to_str(GridWavelength(k) * 1.e3) + ".dat";

      ifstream stream_efficiency_data(file_efficiency_factor.c_str());
      stream_efficiency_data >> description;
      stream_efficiency_data >> description;
      input_tab.Read(stream_efficiency_data, AbsorptionEfficiencyFactor);
      input_tab.Read(stream_efficiency_data, ExtinctionEfficiencyFactor);
      stream_efficiency_data.close();
      cout << "done." << endl;

      for (t = 0; t < Nt_out; t++)
        for (j = 0; j < Ny_out; j++)
          for (i = 0; i < Nx_out; i++)
            {
              for (h = 0; h < Nz_out; h++)
                {
                  real relative_humidity = RelativeHumidity(t, h, j, i);
                  real temperature = Temperature(t, h, j, i);
                  if (relative_humidity > 0.95)
                    relative_humidity = 0.95;
                  if (relative_humidity < 0.05)
                    relative_humidity = 0.05;
                  real coeff_extinction(0.0);
                  real coeff_absorption(0.0);
                  for (is = 0; is < Nbins_in_simulation; is++)
                    {
                      Data<real, 1> index_real(GridSpecies);
                      Data<real, 1> index_imaginary(GridSpecies);
                      real index_dry_real, index_dry_imaginary;
                      real index_wet_real, index_wet_imaginary;
                      Data<real, 1> concentration(GridSpecies);

                      // Black carbon treated as the other components.
                      if (option_black_carbon_treatment == 1)
                        {
                          for (l = 0; l < Nspecies; l++)
                            {
                              concentration(l) = ConcentrationAerosol(l, is,
                                                                      t, h,
                                                                      j, i);
                              index_real(l) = PureSpeciesIndexReal(l, k);
                              index_imaginary(l) =
                                PureSpeciesIndexImaginary(l, k);
                            }

                          if (option_dry_diameter == 1)
                            dry_diameter_computed(is) =
                              compute_Hanel_diameter(relative_humidity,
                                                     wet_diameter);
                          else if (option_dry_diameter == 2)
                            dry_diameter_computed(is) =
                              compute_Gerber_diameter(relative_humidity,
                                                      temperature,
                                                      wet_diameter);
                          else // option_dry_diameter == 3
                            dry_diameter_computed(is) =
                              compute_wet_diameter_from_water_content(
                                                                      concentration.Sum(),
                                                                      ConcentrationWater(is, t, h, j, i),
                                                                      wet_diameter);

                          if (option_wet_index == 1)
                            {
                              if (option_well_mixed_index == 1)
                                // Chemical formula.
                                compute_refractive_index(concentration,
                                                         index_real,
                                                         index_imaginary,
                                                         index_dry_real,
                                                         index_dry_imaginary);
                              else if (option_well_mixed_index == 2)
                                // Lorentz-Lorenz.
                                compute_refractive_index_Lorentz_Lorenz
                                  (concentration,
                                   index_real,
                                   index_imaginary,
                                   index_dry_real,
                                   index_dry_imaginary);

                              compute_Hanel_index(index_dry_real,
                                                  index_dry_imaginary,
                                                  WaterIndexReal(k),
                                                  WaterIndexImaginary(k),
                                                  dry_diameter_computed(is),
                                                  wet_diameter,
                                                  index_wet_real,
                                                  index_wet_imaginary);
                            }
                          else // option_wet_index == 2
                            {
                              Data<real, 1> concentration_tmp(Nspecies + 1);
                              Data<real, 1> index_tmp_real(Nspecies + 1);
                              Data<real, 1> index_tmp_imaginary(Nspecies + 1);
                              for (l = 0; l < Nspecies; l++)
                                {
                                  concentration_tmp(l) = concentration(l);
                                  index_tmp_real(l) = index_real(l);
                                  index_tmp_imaginary(l) = index_imaginary(l);
                                }
                              concentration_tmp(Nspecies) =
                                ConcentrationWater(is, t, h, j, i);
                              index_tmp_real(Nspecies) = WaterIndexReal(k);
                              index_tmp_imaginary(Nspecies) =
                                WaterIndexImaginary(k);

                              if (option_well_mixed_index == 1)
                                // Chemical formula.
                                compute_refractive_index(concentration_tmp,
                                                         index_tmp_real,
                                                         index_tmp_imaginary,
                                                         index_wet_real,
                                                         index_wet_imaginary);
                              else if (option_well_mixed_index == 2)
                                // Lorentz-Lorenz.
                                compute_refractive_index_Lorentz_Lorenz
                                  (concentration_tmp,
                                   index_tmp_real,
                                   index_tmp_imaginary,
                                   index_wet_real,
                                   index_wet_imaginary);
                            }
                        }
                      else if (option_black_carbon_treatment == 2)
                        // Black carbon core treatement.
                        {
                          real index_wet_no_black_carbon_real;
                          real index_wet_no_black_carbon_imaginary;
                          Data<real, 1> concentration_no_black_carbon
                            (GridSpeciesNoBlackCarbon);
                          Data<real, 1> concentration_water_no_black_carbon
                            (GridSpecies);
                          real concentration_black_carbon =
                            ConcentrationAerosol(black_carbon_index, is, t,
                                                 h, j, i);
                          int ll = 0;
                          for (l = 0; l < Nspecies; l++)
                            {
                              concentration(l) =
                                ConcentrationAerosol(l, is, t, h, j, i);
                              if (l != black_carbon_index)
                                {
                                  concentration_no_black_carbon(ll) =
                                    ConcentrationAerosol(l, is, t, h, j, i);
                                  concentration_water_no_black_carbon(ll) =
                                    ConcentrationAerosol(l, is, t, h, j, i);
                                  index_real(ll) = PureSpeciesIndexReal(l, k);
                                  index_imaginary(ll) =
                                    PureSpeciesIndexImaginary(l, k);
                                  ll += 1;
                                }
                            }

                          concentration_water_no_black_carbon(Nspecies - 1) =
                            ConcentrationWater(is, t, h, j, i);

                          index_real(Nspecies - 1) = WaterIndexReal(k);
                          index_imaginary(Nspecies - 1) =
                            WaterIndexImaginary(k);

                          if (ConcentrationWater(is, t, h, j, i) < 1.e4)
                            {
                              compute_wet_diameter_from_water_content
                                (dry_diameter_computed(is),
                                 concentration_no_black_carbon.Sum(),
                                 ConcentrationWater(is, t, h, j, i),
                                 wet_diameter_without_black_carbon);
                              compute_wet_diameter_from_water_content
                                (dry_diameter_computed(is),
                                 concentration.Sum(),
                                 ConcentrationWater(is, t, h, j, i),
                                 wet_diameter);
                            }
                          else
                            {
                              compute_Hanel_diameter
                                (dry_diameter_computed(is),
                                 relative_humidity,
                                 wet_diameter);
                              wet_diameter_without_black_carbon
                                = wet_diameter;
                            }
                          // Solution without black carbon.
                          if (option_well_mixed_index == 1)
                            // Chemical formula.
                            compute_refractive_index
                              (concentration_water_no_black_carbon,
                               index_real,
                               index_imaginary,
                               index_wet_no_black_carbon_real,
                               index_wet_no_black_carbon_imaginary);
                          else if (option_well_mixed_index == 2)
                            // Lorentz-Lorenz.
                            compute_refractive_index_Lorentz_Lorenz
                              (concentration_water_no_black_carbon,
                               index_real,
                               index_imaginary,
                               index_wet_no_black_carbon_real,
                               index_wet_no_black_carbon_imaginary);

                          compute_refractive_index_Maxwell_Garnet
                            (concentration_black_carbon,
                             concentration_water_no_black_carbon.Sum(),
                             PureSpeciesIndexReal(black_carbon_index, k),
                             PureSpeciesIndexImaginary(black_carbon_index, k),
                             index_wet_no_black_carbon_real,
                             index_wet_no_black_carbon_imaginary,
                             index_wet_real, index_wet_imaginary);
                        }
                      RegularGrid<real> grid_tmp_index_real(1);
                      RegularGrid<real> grid_tmp_index_imaginary(1);
                      RegularGrid<real> grid_tmp_diameter(1);
                      grid_tmp_index_real(0) = index_wet_real;
                      grid_tmp_index_imaginary(0) = index_wet_imaginary;
                      grid_tmp_diameter(0) = wet_diameter;
                      Data<real, 3> data_absorption_efficiency_factor
                        (grid_tmp_index_real, grid_tmp_index_imaginary,
                         grid_tmp_diameter);
                      Data<real, 3> data_extinction_efficiency_factor
                        (grid_tmp_index_real, grid_tmp_index_imaginary,
                         grid_tmp_diameter);
                      real efficiency_absorption_factor;

                      LinearInterpolationRegular
                        (AbsorptionEfficiencyFactor,
                         data_absorption_efficiency_factor);

                      LinearInterpolationRegular
                        (ExtinctionEfficiencyFactor,
                         data_extinction_efficiency_factor);

                      efficiency_absorption_factor =
                        (data_extinction_efficiency_factor(0, 0, 0) -
                         data_absorption_efficiency_factor(0, 0, 0));

                      if (concentration.Sum() > 0.0)
                        {
                          coeff_extinction += 3.
                            * data_extinction_efficiency_factor(0, 0, 0)
                            * concentration.Sum()
                            * wet_diameter * wet_diameter
                            / (2. * dry_aerosol_density
                               * dry_diameter_computed(is)
                               * dry_diameter_computed(is)
                               * dry_diameter_computed(is));

                          coeff_absorption += 3.
                            * efficiency_absorption_factor
                            * concentration.Sum()
                            * wet_diameter * wet_diameter
                            / (2. * dry_aerosol_density
                               * dry_diameter_computed(is)
                               * dry_diameter_computed(is)
                               * dry_diameter_computed(is));
                        }
                    }

                  ExtinctionCoefficient(k, t, h, j, i) = coeff_extinction;
                  AbsorptionCoefficient(k, t, h, j, i) = coeff_absorption;
                  OpticalThickness(k, t, j, i) += coeff_extinction
                    * (GridZ_interf_out(h + 1) - GridZ_interf_out(h));
                  OpticalThicknessScattering(k, t, j, i) +=
                    (coeff_extinction - coeff_absorption)
                    * (GridZ_interf_out(h + 1) - GridZ_interf_out(h));
                }

              SingleScatteringAlbedo(k, t, j, i) =
                OpticalThicknessScattering(k, t, j, i)
                / OpticalThickness(k, t, j, i);
            }
      cout << "done." << endl;
    }

  cout << "Writing output results...";
  cout.flush();

  FormatBinary<float> output;
  for (k = 0; k < Nwavelength; k++)
    {
      Array<real, 3> OpticalThickness_extract = OpticalThickness.GetArray()
        (k, Range::all(), Range::all(), Range::all());

      output.Append(OpticalThickness_extract,
                    directory_result + string("OpticalThickness_")
                    + to_str(GridWavelength(k) * 1.e3) + ".bin");
    }

  for (k = 0; k < Nwavelength; k++)
    {
      Array<real, 4> coeffext_extract = ExtinctionCoefficient.GetArray()
        (k, Range::all(), Range::all(), Range::all(), Range::all());

      output.Append(coeffext_extract, directory_result
                    + string("ExtinctionCoefficient_")
                    + to_str(GridWavelength(k) * 1.e3) + ".bin");

    }

  for (k = 0; k < Nwavelength; k++)
    {
      Array<real, 4> coeffabs_extract = AbsorptionCoefficient.GetArray()
        (k, Range::all(), Range::all(), Range::all(), Range::all());

      output.Append(coeffabs_extract, directory_result
                    + string("AbsorptionCoefficient_")
                    + to_str(GridWavelength(k) * 1.e3) + ".bin");
    }

  for (k = 0; k < Nwavelength; k++)
    {
      Array<real, 3> ssa_extract = SingleScatteringAlbedo.GetArray()
        (k, Range::all(), Range::all(), Range::all());

      output.Append(ssa_extract, directory_result
                    + string("SingleScatteringAlbedo_")
                    + to_str(GridWavelength(k) * 1.e3) + ".bin");
    }
  cout << " done." << endl;

  cout << endl;

  END;

  return 0;

}
