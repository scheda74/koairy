// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Yelva Roustan
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


#ifndef POLYPHEMUS_FILE_MODELS_BASEPOINTEMISSION_HXX



//////////////
// INCLUDES //
//////////////


#include "ContinuousEmission.hxx"
#include "ContinuousEmission_aer.hxx"
#include "ContinuousLineEmission.hxx"
#include "ContinuousLineEmission_aer.hxx"
#include "PuffEmission.hxx"
#include "PuffEmission_aer.hxx"
#include "TemporalEmission.hxx"
#include "TemporalEmission_aer.hxx"


namespace Polyphemus
{
  using namespace std;
  using namespace AtmoData;

  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages a set of point emissions.
  template<class T>
  class BasePointEmission
  {

  protected:

    //! List of point emissions.
    vector<PointEmissionUnit<T>*> Emission;
    //! Source id
    string source_id;

  public:

    BasePointEmission();
    virtual ~BasePointEmission();

    virtual void Init(string config_file, vector<string> species_list);
    virtual void Init(string config_file, vector<string> species_list,
                      vector<string> species_list_aer, int Nbin_aer);

    virtual int GetNumberEmission();
    string GetEmissionSourceId(int emission);
    T GetEmissionSourceWater(int emission) const;
    void GetEmissionCoordinates(T& abscissa, T& ordinate, T& height,
                                int emission) const;
    void GetEmissionCoordinates(list<Array<T, 1> >& coordinate_list,
                                int emission);
    void SetEmissionCoordinates(T abscissa, T ordinate, T height,
                                int emission);
    void SetEmissionCoordinates(const list<Array<T, 1> >& coordinate_list,
                                int emission);
    string GetEmissionType(int emission);
    Date GetDateBeg(int emission);
    Date GetDateEnd(int emission);
    bool HasBegun(Date date, int emission);
    bool HasEnded(Date date, int emission);
    vector<int> GetEmittedSpeciesIndex(int emission) const;
    vector<int> GetEmittedSpeciesIndex_aer(int emission) const;
    vector<map <string, string> > GetEmittedAerosolSpecies(int emission) const;

    T GetRate(int emission, int species) const;
    T GetRate(int emission, int species, int id_section) const;
    void MultiplyRate(T factor);
    void MultiplyRate(int species, T factor);
    void ApplyTimeShift(T shift);
    void ApplyAltitudeShift(T shift);
    void GetEmission(Date date_beg, Date date_end, int species,
                     int emission, Array<T, 2>& point_emission);
    void GetEmission_aer(Date date_beg, Date date_end,
                         int species, int emission,
                         T& quantity_aer_bin);
    bool IsEmitting(Date current_date, Date next_date, int emission);
    bool HasPlumeRise(int emission);
    void GetPlumeRiseParam(T& velocity, T& temperature,
                           T& diameter, int emission);
    T GetEmissionTimespan(const Date& current_date,
                          const Date& next_date, int emission) const;
    T GetWidth(int emission, int id_section);
    T GetVehicleVelocity(int emission, int id_section);
    T GetArea(int emission, int id_section);
    T GetDensity(int emission, int id_section);
    void AddContinuousEmission(T abscissa, T ordinate, T height,
                               const vector<T>&  rate, Date date_beg,
                               Date date_end, const vector<int>& species,
                               T diameter, T velocity, T temperature);

    bool IsVolumeSource(int emission);
    void GetVolumeSource(T& width, T& length, int emission);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_BASEPOINTEMISSION_HXX
#endif
