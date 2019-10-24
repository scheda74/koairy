// Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
// Author(s): Meryem Ahmed de Biasi
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


#ifndef POLYPHEMUS_FILE_MODELS_STATIONARYMODEL_HXX


#include "BaseModel.cxx"
#include <map>
#include <vector>
#include <string>
#include "AtmoData.hxx"


namespace Polyphemus
{


  using namespace std;
  using namespace AtmoData;


  //////////////////////
  // STATIONARYMODEL //
  ////////////////////


  //! Model that performs a time-integration at local scale.
  /*!
    This model is mainly used as an interface between the driver and the
    underlying model which is used to perform an inner loop to reach
    convergence.
  */
  template<class T, class ClassModel>
  class StationaryModel: public BaseModel<T>
  {

  protected:

    //! Underlying model.
    ClassModel Model;

    /*** Convergence ***/

    //! Convergence criterion.
    T epsilon;
    //! Norm used to check convergence.
    string norm;
    //! Method used to normalize norm (with the mean or the maximum value).
    string method;
    //! Bool which evaluates if convergence has been reached.
    bool convergence;

    // Data used to check if convergence has been reached.
    //! Data storing concentrations before the iteration.
    Data<T, 4> Conc_prev;
    //! Data storing concentrations after the iteration.
    Data<T, 4> Conc_current;

  public:

    /*** Constructors and destructor ***/

    StationaryModel(string config_file);
    virtual ~StationaryModel();

    /*** Other methods ***/

    void Init();
    void InitStep();
    void Forward();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STATIONARYMODEL_HXX
#endif
