// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
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

// This file is part of the Eulerian model Castor.

// This code is essentially based on the chemistry-transport model Chimere,
// distributed under GNU GPL -- copyright (C) 2005 Institut Pierre-Simon
// Laplace (CNRS), INERIS, LISA (CNRS).


#ifndef POLYPHEMUS_FILE_MODELS_CASTORCHEMISTRY_HXX


#include <vector>
#include "AtmoData.hxx"
#include "CastorTransport.cxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  /////////////////////
  // CASTORCHEMISTRY //
  /////////////////////


  //! This class is a solver for an advection-diffusion-reaction equation.
  /*!  CastorChemistry is a clone of Chimere with respect to transport and
    gaseous chemistry.
  */
  template<class T, class ClassTransport, class ClassChemistry>
  class CastorChemistry:
    public CastorTransport<T, ClassTransport>
  {

    /*** Output ***/

    //! Output unit.
    string unit;
    //! Molecular weight of species (g / mol).
    map<string, T> molecular_weight;

  public:

    /*** Type declarations ***/

    typedef typename map<string, InputFiles<T> >::iterator
    input_files_iterator;

  public:

    /*** Humidity ***/

    //! Specific humidity at current date.
    Data<T, 3> SpecificHumidity;
    //! Specific humidity buffer.
    Data<T, 3> FileSpecificHumidity_i;
    //! Specific humidity buffer.
    Data<T, 3> FileSpecificHumidity_f;

    /*** Attenuation ***/

    //! Attenuation coefficients at current date.
    Data<T, 2> Attenuation;
    //! Attenuation buffer.
    Data<T, 2> FileAttenuation_i;
    //! Attenuation buffer.
    Data<T, 2> FileAttenuation_f;

    //! Liquid water content (kg / kg).
    Data<T, 3> LiquidWaterContent;
    //! Liquid water content buffer.
    Data<T, 3> FileLiquidWaterContent_i;
    //! Liquid water content buffer.
    Data<T, 3> FileLiquidWaterContent_f;

    /*** Chemical mechanism ***/

    //! Chemical mechanism and chemical numerical scheme.
    ClassChemistry Chemistry_;

  public:

    /*** Constructor and destructor ***/

    CastorChemistry(string config_file);
    virtual ~CastorChemistry();

    /*** Configuration ***/

    virtual void ReadConfiguration();
    virtual void CheckConfiguration();

    int GetSpeciesIndex_ext(string species) const;

    /*** Initializations ***/

    virtual void Allocate();
    void Init();
    void InitStep();

    virtual void SetDate(Date date);

    /*** Integration ***/

    void Forward();

    /*** Units conversion ***/

    void ConvertExternalUnit(Data<T, 4>& Concentration_);
    void ConvertInternalUnit(Data<T, 4>& Concentration_);

  protected:

    virtual void InitAllData();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CASTORCHEMISTRY_HXX
#endif
