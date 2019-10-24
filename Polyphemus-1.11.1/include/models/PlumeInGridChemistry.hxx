// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok
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


#ifndef POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDCHEMISTRY_HXX

#include "PlumeInGridTransport.cxx"

namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////////
  // PLUME-IN-GRID DRIVER //
  //////////////////////////


  /*! \brief This class provides a driver to perform a plume-in-grid
    simulation.
  */
  /*! The driver is responsible for the model initialization, for the loop
    over all meteorological conditions and over all time steps, and for the
    calls to the output saver. It also Its reference floating-point precision
    is 'T'. The model is an instance of 'ClassModel' and the output saver is
    an instance of 'ClassOutputSaver'.
  */
  template<class T, class ClassEulerianModel, class ClassLocalModel>
  class PlumeInGridChemistry:
    public PlumeInGridTransport<T, ClassEulerianModel, ClassLocalModel>
  {

  protected:

    /*** Main components ***/

    //! Is there chemistry?
    bool option_chemistry;


    //! Is there feedback for chemistry?
    bool option_chemical_feedback;

    //! Number of species with photolysis.
    int Nr_photolysis;

    /*** Meteorological data for chemistry ***/

    //! Attenuation.
    Data<T, 3> Attenuation_i;
    //! Pressure.
    Data<T, 3> Pressure_i;
    //! Specific humidity.
    Data<T, 3> SpecificHumidity_i;

    /*** Other data for chemistry processes ***/

    //! Photolysis rates.
    Data<T, 4> PhotolysisRate_i;

  public:

    /*** Constructor and destructor ***/

    PlumeInGridChemistry(string config_file);
    virtual ~PlumeInGridChemistry();

    /*** General methods ***/

    void ReadConfiguration();
    void CheckConfiguration();
    void Init();
    void InitStep();
    void Forward();

    /*** Meteorological and species data ***/

    virtual void UpdateMeteo(int puff_index);
    void ExtractMeteo(T height, T lat, T lon, bool isday,
                      bool option_similarity, map<string, T>& met_data,
                      string& stability, bool& rural);
    void ExtractSpeciesData(T height, T lat, T lon,
                            Array<T, 1>& deposition_velocity,
                            Array<T, 1>& scavenging_coefficient,
                            Array<T, 1>& photolysis_rate,
                            Array<T, 1>& background_concentration);

    /*** Puff transfer ***/

    void ComputeChemistry(Array<int, 2> PuffCellList,
                          Array<T, 4>& PuffConcentration);

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_DRIVER_PLUMEINGRIDCHEMISTRY_HXX
#endif
