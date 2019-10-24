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
// For more information, visit the Polyphemus web site:
//      http://cerea.enpc.fr/polyphemus/


#ifndef POLYPHEMUS_FILE_DRIVER_STREETINGRIDCHEMISTRY_HXX

#include "StreetInGridTransport.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  ////////////////////////////////////
  // STREET-IN-GRID CHEMISTRY MODEL //
  ////////////////////////////////////


  /*! \brief This class provides an interface to perform a street-in-grid
    simulation.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  class StreetInGridChemistry:
    public StreetInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  {

  protected:

    /*** Main components ***/

    //! Is there chemistry?
    bool option_chemistry;

    string chemical_mechanism;

    string chemical_mechanism_local;

    //! Number of species with photolysis.
    int Nr_photolysis;

    /*** Meteorological data for chemistry ***/

    //! Attenuation.
    Data<T, 3> Attenuation_i;
    //! Pressure.
    Data<T, 3> Pressure_i;
    //! Specific humidity.
    Data<T, 3> SpecificHumidity_i;
    //! Temperature.
    Data<T, 3> Temperature_i;

    /*** Other data for chemistry processes ***/

    //! Photolysis rates.
    Data<T, 4> PhotolysisRate_i;

  public:

    /*** Constructor and destructor ***/

    StreetInGridChemistry(string config_file);
    virtual ~StreetInGridChemistry();

    /*** General methods ***/

    void ReadConfiguration();
    void CheckConfiguration();
    void Allocate();
    void Init();
    void InitStep();
    void Forward();

    /*** Meteorological and species data ***/

    void UpdateStreetData(int street_index);
    void ExtractAdditionalMeteo(T height, T lat, T lon, map<string, T>& met_data);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_STREETINGRIDCHEMISTRY_HXX
#endif
