// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Ir�ne Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_AER_HXX



//////////////
// INCLUDES //
//////////////


#include "PointEmissionUnit.cxx"

namespace Polyphemus
{
  using namespace std;
  using namespace AtmoData;

  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages point emissions.
  template<class T>
  class PointEmissionUnit_aer: public PointEmissionUnit<T>
  {

  protected:

    // //! Type of emission.
    // string type;



    // //! List of emitted species.
    // vector<int> emitted_species_index;
    vector<int> emitted_species_index_aer;
    vector<map <string, string> > emitted_species_list_aer_bin;

    // //! Number of emitted species.
    // int Ns_emis;
    int Ns_emis_aer;

    // //! Beginning date of emission.
    // Date date_beg;

    // //! Efflux speed of gases.
    // T velocity_;

    // //! Temperature.
    // T temperature_;

    // //! Diameter.
    // T diameter_;


  public:

    PointEmissionUnit_aer();
    virtual ~PointEmissionUnit_aer();
    virtual void Init(ConfigStream& config, vector<string> species_list, vector<string> species_list_aer);
    // void SetType(string type_name);
    // virtual void SetCoordinates(T abscissa, T ordinate, T height);
    // string GetType() const;
    // virtual void GetCoordinates(T& abscissa, T& ordinate, T& height) const;
    // virtual void GetCoordinates(list<Array<T, 1> >& coordinate_list) const;
    // Date GetDateBeg() const;
    // vector<int> GetEmittedSpeciesIndex() const;
    int GetSpeciesIndex_aer(string species,
                            vector<string> species_list_aer) const;
    // virtual T GetEmissionTimespan(const Date& current_date,
    //            const Date& next_date) const;
    // bool HasBegun(Date date);
    // virtual bool HasEnded(Date date);
    // virtual T GetRate(int species) const;
    // virtual void GetEmission(Date date_beg, Date date_end, int s,
    //             Array<T, 2>& point_emission_list);
    // virtual void GetPlumeRiseParam(T& velocity, T& temperature,
    //             T& diameter);
    // bool HasPlumeRise();
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_AER_HXX
#endif
