// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Ir�ne Korsakissok, Yelva Roustan, Youngseob Kim
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


#ifndef POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_HXX



//////////////
// INCLUDES //
//////////////


#include <vector>
#include <map>
#include <fstream>

#include "AtmoDataHeader.hxx"


namespace Polyphemus
{
  using namespace std;
  using namespace AtmoData;

  /////////////////////
  // POINT EMISSIONS //
  /////////////////////


  //! This class manages point emissions.
  template<class T>
  class PointEmissionUnit
  {

  protected:

    //! Type of emission.
    string type;

    //! List of emitted species.
    vector<int> emitted_species_index;
    vector<int> emitted_species_index_aer;
    vector<map <string, string> > emitted_species_list_aer_bin;

    //! Number of emitted species.
    int Ns_emis;
    int Ns_emis_aer;

    //! Beginning date of emission.
    Date date_beg;
    Date date_end;

    //! Efflux speed of gases.
    T velocity_;

    //! Temperature.
    T temperature_;

    //! Diameter.
    T diameter_;

    //! Emitted Water.
    T source_water;

    //! Source id
    string source_id;

    //! if the source is a volume source.
    bool is_volume_source_;
    T width_;
    T length_;

  public:

    PointEmissionUnit();
    virtual ~PointEmissionUnit();

    virtual void Init(ConfigStream& config, vector<string> species_list);
    virtual void Init(ConfigStream& config, vector<string> species_list, vector<string> species_list_aer, int Nbin_aer);
    virtual void Init(T abscissa, T ordinate, T height,
                      const vector<T>& rate, Date date_beg, Date date_end,
                      const vector<int>& species, T diameter, T velocity,
                      T temperature);

    void SetType(string type_name);
    string GetSourceId() const;
    T GetSourceWater() const;
    virtual void SetCoordinates(T abscissa, T ordinate, T height);
    virtual void SetCoordinates(const list<Array<T, 1> >& coordinate_list);
    string GetType() const;
    virtual void GetCoordinates(T& abscissa, T& ordinate, T& height) const;
    virtual void GetCoordinates(list<Array<T, 1> >& coordinate_list) const;
    Date GetDateBeg() const;
    Date GetDateEnd() const;
    vector<int> GetEmittedSpeciesIndex() const;
    vector<int> GetEmittedSpeciesIndex_aer() const;
    vector<map <string, string> > GetEmittedAerosolSpecies() const;

    int GetSpeciesIndex(string species,
                        vector<string>& ref_species_list) const;
    int GetSpeciesIndex_aer(string species,
                            vector<string> species_list_aer) const;

    virtual T GetEmissionTimespan(const Date& current_date,
                                  const Date& next_date) const;
    bool HasBegun(Date date);
    virtual bool HasEnded(Date date);
    virtual T GetRate(int species) const;
    virtual T GetRate(int species, int id_section) const;
    virtual void MultiplyRate(T factor);
    virtual void MultiplyRate(int species, T factor);
    virtual void ApplyTimeShift(T shift);
    virtual void ApplyAltitudeShift(T shift);
    virtual void GetEmission(Date date_beg, Date date_end, int s,
                             Array<T, 2>& point_emission_list);
    virtual void GetEmission_aer(Date date_beg, Date date_end, int s, T& quantity);
    virtual void GetPlumeRiseParam(T& velocity, T& temperature,
                                   T& diameter);

    bool HasPlumeRise();
    virtual T GetWidth(int id_section) const;
    virtual int GetIdSource() const;
    virtual int GetIdSection() const;
    virtual T GetVehicleVelocity(int id_section) const;
    virtual T GetArea(int id_section) const;
    virtual T GetDensity(int id_section) const;

    bool IsVolumeSource();
    virtual void GetVolumeSource(T& width, T& length);
  };

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_POINTEMISSIONUNIT_HXX
#endif
