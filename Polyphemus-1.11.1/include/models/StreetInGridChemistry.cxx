// Copyright (C) 2016, CEREA - ENPC - EDF R&D
// Author(s): Youngseob Kim
//
// This file is part of the air quality modeling system Polyphemus.
//
// Polyphemus is developed in the ENPC - EDF R&D joint laboratory CEREA.
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
// For more information, visit the Polyphemus web site
//      http://cerea.enpc.fr/polyphemus/


#ifndef POLYPHEMUS_FILE_MODELS_STREETINGRIDCHEMISTRY_CXX


#include "StreetInGridChemistry.hxx"


namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////


  //! Main constructor.
  /*! Builds the model and reads option keys in the configuration file.
    \param config_file configuration file.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::StreetInGridChemistry(string config_file): 
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>(config_file)
  {
  }


  //! Destructor.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::~StreetInGridChemistry()
  {
  }


  /////////////////////
  // INITIALIZATIONS //
  /////////////////////


  //! Reads the configuration.
  /*! It reads the description of the domain, the simulation starting-date,
    the species list and display options.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ReadConfiguration()
  {
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
    ::ReadConfiguration();

    this->config.SetSection("[options]");
    bool option_photolysis;
    this->config.PeekValue("With_chemistry", option_chemistry);
    this->config.PeekValue("With_photolysis", option_photolysis);
    option_chemistry = (option_chemistry && option_photolysis);
    
    this->config.PeekValue("Option_chemistry", chemical_mechanism);
  }

  //! Checks that the configuration is acceptable.
  /*! In case any inconsistency is found, an exception is thrown.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::CheckConfiguration()
  {
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::CheckConfiguration();
  }

  //! Allocates memory.
  /*! Allocates grids and fields.
   */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::Allocate()
  {
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
      ::Allocate();
  }


  //! Initialization.
  /*! Initializes the Eulerian model and the meteorological conditions. Reads
    the sources and creates the corresponding Gaussian model instances.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::Init()
  {
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
    ::Init();

    /*** Initializations ***/

    Nr_photolysis = 0;
    vector<string> photolysis_reaction_list;

    if (option_chemistry)
      {
	// Extracts information useful for chemistry from Eulerian model.
	Nr_photolysis = this->Model.GetNr_photolysis();
	photolysis_reaction_list = this->Model.GetPhotolysisReactionList();
	
	// Photolysis rates.
	if (this->Model.HasField("PhotolysisRate"))
	  PhotolysisRate_i.Copy(this->Model.D4("PhotolysisRate_i"));
	else
	  throw string("Eulerian model doesn't have a field ")
	    + "\"PhotolysisRate\".";

	// Meteorological data.
	Attenuation_i.Copy(this->Model.D3("Attenuation_i"));
	Pressure_i.Copy(this->Model.D3("Pressure_i"));
	SpecificHumidity_i.Copy(this->Model.D3("SpecificHumidity_i"));
	Temperature_i.Copy(this->Model.D3("Temperature_i"));
      }

    // Get local chemical mechanism.
    if (option_chemistry)
      {
        chemical_mechanism_local = this->StreetNetworkModel->GetChemicalMechanism();
        if (chemical_mechanism != chemical_mechanism_local)
          {
            cout << "Warning: chemical mechanisms in the host model and the local model are different: host (" << chemical_mechanism << ") and local (" << chemical_mechanism_local << ")" << endl;
          }
      }


  }


  //! Model initialization for each step.
  /*! It reads on file the data that are is needed for the current step.
   */  
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::InitStep()
  {
    //! InitStep of Polair3DChemistry
    this->Model.InitStep();

    //! UpdataMeteo of StreetInGridTransport
    this->UpdateMeteo();

    if (option_chemistry)
      {
        Attenuation_i.Copy(this->Model.D3("Attenuation_i"));
        Pressure_i.Copy(this->Model.D3("Pressure_i"));
        SpecificHumidity_i.Copy(this->Model.D3("SpecificHumidity_i"));
        Temperature_i.Copy(this->Model.D3("Temperature_i"));

        // Photolysis rates.
        PhotolysisRate_i.Copy(this->Model.D4("PhotolysisRate_i"));
      }


    /*** Meteo data on the streets***/

    // Updating the meteorological data for the streets.
    for (int street_index = 0; street_index < this->Nstreet; street_index++)
      UpdateStreetData(street_index);

    // Updating the meteorological data for the intersections.
    this->Nintersection = this->StreetNetworkModel->GetNumberIntersection();
    for (int intersection_index = 0; intersection_index < this->Nintersection;
         intersection_index++)
      this->UpdateIntersectionData(intersection_index);

    this->StreetNetworkModel->InitStep();

    bool emission_save = false;
    if (emission_save)
      this->EmissionRateSave();
  }


  //! Performs one step forward.
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::Forward()
  {
    int rank;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = 0;
#endif

    if (rank == 0)
      {
        this->cell_quantity.resize(this->Ns, this->Ny, this->Nx);
        this->cell_quantity = 0.0;

        this->CorrectBackgroundConcentration();

        this->UpdateBackgroundConcentration();

        // Street-network model.
        this->StreetNetworkModel->Transport();

        if (option_chemistry)
          this->StreetNetworkModel->Chemistry();

        this->StreetNetworkModel->SetStreetConcentration();

        // Compute the mass flux.
        this->ComputeMassFlux();

        // Correct the background concentration by the ratio of the urban volume.
        this->RestoreBackgroundConcentration();
      }

    // Eulerian model.
    this->Model.SetForward();

    this->Model.AdvectionSplit();

    this->Model.DiffusionSplitXY();

    if (rank == 0)
      this->AddMassFluxToBackgroundConcentration();

    this->Model.DiffusionSplitZ();

    this->Model.ChemistrySplit();
 
    if (rank == 0)
      {

        /* Adding all streets to concentration Data in case concentrations are saved
           on the domain.*/

        this->StreetConcentration.Copy(this->StreetNetworkModel->GetStreetConcentration());

        this->Concentration.Copy(this->Model.GetConcentration());
        this->Nstreet = this->StreetNetworkModel->GetNStreet();
        for (int street_index = 0; street_index < this->Nstreet; street_index++)
          {
            // Street center coordinates.
            T z_c, lon_c, lat_c;
            z_c = 2.0;
            this->StreetNetworkModel->
              GetStreetCoordinate(street_index, lon_c, lat_c);
            int iz, iy, ix;
            this->GetCellIndices(lon_c, lat_c, z_c, iz, iy, ix);

            for (int s = 0; s < this->Ns; ++s)
              {
                T quantity_street = this->StreetNetworkModel->
                  GetStreetQuantity(street_index, s);
                this->cell_quantity(s, iy, ix) += quantity_street;
              }
          }
	
        //! Adding street concentration to the Concentration Data.
        this->StreetTransfer(this->Concentration);

        /* Increase the time step */
        this->StreetNetworkModel->AddTime(this->delta_t_local);
        this->AddTime(this->delta_t_eulerian);
        this->step++;
      }
  }

  //! Updates the street meteorological data.
  /*! Updates the street meteorological data taking them from the center cell.
    \param street_index index of the street.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::UpdateStreetData(int street_index)
  {
    StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
    ::UpdateStreetData(street_index);

     if (option_chemistry)
       { 
         // Street center coordinates.
         T lon_c, lat_c;
         T z_c = 2.0; 
         this->StreetNetworkModel->
           GetStreetCoordinate(street_index, lon_c, lat_c);
    
         // Extracting meteorological data.
         map<string, T> meteorological_data;
         ExtractAdditionalMeteo(z_c, lat_c, lon_c, meteorological_data);

         this->StreetNetworkModel->
           SetStreetAdditionalMeteo(street_index,
                                    meteorological_data["attenuation"],
                                    meteorological_data["pressure"],
                                    meteorological_data["specific_humidity"],
                                    meteorological_data["temperature"]);
       }
  }

  //! Extracts meteorological data.
  /*! Extracts meteorological data from 3D fields at a given point.
    \param D3_map pointer to map containing 3D meteorological data.
    \param index_z cell index along z.
    \param index_y cell index along y.
    \param index_x cell index along x.
    \param met_data map containing meteorological data at the given point.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  void StreetInGridChemistry<T, ClassEulerianModel, ClassLocalModel>
  ::ExtractAdditionalMeteo(T height, T lat, T lon, map<string, T> &met_data)
  {
    int index_x, index_y, index_z;
    this->GetCellIndices(lon, lat, height, index_z, index_y, index_x);
    Array<T, 1> Coord3D(3);
    Coord3D(0) = height;
    Coord3D(1) = lat;
    Coord3D(2) = lon;

    T attenuation, pressure, specific_humidity, temperature;

    // Test if other fields have to be interpolated.
    if (!this->interpolated)
      {
	attenuation = Attenuation_i(index_z, index_y, index_x);
	pressure = Pressure_i(index_z, index_y, index_x);
	specific_humidity = SpecificHumidity_i(index_z, index_y, index_x);
	temperature = Temperature_i(index_z, index_y, index_x);
      }
    else
      {
        LinearInterpolationPoint(Attenuation_i, Coord3D, attenuation);
        LinearInterpolationPoint(Pressure_i, Coord3D, pressure);
        LinearInterpolationPoint(SpecificHumidity_i, Coord3D, specific_humidity);
        LinearInterpolationPoint(Temperature_i, Coord3D, temperature);
      }

    met_data["attenuation"] = attenuation;
    met_data["pressure"] = pressure;
    met_data["specific_humidity"] = specific_humidity;
    met_data["temperature"] = temperature;
  }
  
} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETINGRIDCHEMISTRY_CXX
#endif
