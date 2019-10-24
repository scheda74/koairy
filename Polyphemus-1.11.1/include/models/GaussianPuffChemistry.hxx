// Copyright (C) 2005-2012, ENPC - INRIA - EDF R&D
// Author(s): Irène Korsakissok, Vivien Mallet, Youngseob Kim
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

// This file is part of a Gaussian puff model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFCHEMISTRY_HXX


#include "PuffChemistry.cxx"
#include "GaussianPuffTransport.cxx"


namespace Polyphemus
{


  using namespace std;


  /////////////////////////////
  // GAUSSIAN PUFF CHEMISTRY //
  /////////////////////////////


  //! Gaussian puff model with gaseous chemistry.
  template<class T, class ClassChemistry>
  class GaussianPuffChemistry: public GaussianPuffTransport<T>
  {

  protected:


    /********* Chemistry *********/

    //! Number of subcycles for chemistry.
    int Ncycle;

    //! File to save the total plume mass.
    string file_mass;

    //! List of photolysis rates.
    Array<T, 1> photolysis_rate_;

    //! List of background concentrations.
    Array<T, 1> background_concentration_;

    //! List of species with photolysis reactions.
    vector<string> photolysis_reaction_list;

    //! Number of photolysis reactions.
    int Nr_photolysis;

    //! List of species with forced concentrations.
    vector<string> species_list_forced;

    //! Chemical mechanism and chemical numerical scheme.
    ClassChemistry Chemistry_;

    /********* Background Meteorological data ************/

    //! Attenuation.
    T background_attenuation_;
    //! Pressure.
    T background_pressure_;
    //! Specific humidity.
    T background_specific_humidity_;

    //! longitude;
    T background_longitude;
    //! latitude;
    T background_latitude;

    /********* Current Meteorological data ************/

    //! Attenuation.
    T attenuation_;
    //! Pressure.
    T pressure_;
    //! Specific humidity.
    T specific_humidity_;
    //! Specific humidity tmp.
    T specific_humidity_tmp;
    //! Puff Water.
    T water_puff_;

  public:

    /*** Constructors and destructor ***/

    GaussianPuffChemistry(string config_file);
    virtual ~GaussianPuffChemistry();

    /*** Initializations ***/

    virtual void ReadConfiguration();
    virtual void Allocate();
    void Init();
    void InitStep();
    virtual void InitPuff();

    /*** Access Methods ***/

    int GetNr_photolysis() const
    {
      return Nr_photolysis;
    }
    int GetNs_source() const
    {
      return 0;
    }
    int GetNz_source() const
    {
      return 0;
    }
    int SourceGlobalIndex(int s) const
    {
      return 0;
    }
    vector<string> GetPhotolysisReactionList() const
    {
      return photolysis_reaction_list;
    }
    bool WithChemistry();
    bool WithPhotolysis();

    /*** Background Meteorological data ***/

    virtual void InitMeteo(ConfigStream& meteo, bool show_meteo);
    void InitPhotolysis(int Nr, vector<string> photolysis_list);
    virtual void InitChemistry(ConfigStream& meteo);
    virtual void SetCurrentMeteo(Puff<T>* puff);

    /*** Access methods for Puff attributes ***/

    void SetPuffCurrentMeteo(T interaction_coefficient, Puff<T>* puff);
    void SetPuffAdditionalMeteo(int index, T attenuation, T pressure,
                                T specific_humidity);
    void SetPuffPhotolysisRate(int index, Array<T, 1> photolysis_rate);
    void SetPuffBackgroundConcentration(int index, Array<T, 1> concentration);

    T GetPuffBackgroundConcentration(int index, int s);
    void SetPuffQuantity(int index, Array<T,1> quantity);
    T GetPuffQuantity(int index, int s);


    T GetSpecificHumidity();


    /*** Computational Methods ***/

    T GetConcentration(int species, T z, T y, T x);
    T GetIntegratedConcentration(int species, T z, T y, T x,
                                 T lz, T ly, T lx);
    using GaussianPuffTransport<T>::GetConcentration;

    virtual void Chemistry();
    virtual void Chemistry(vector<list<int> > PuffCellList,
                           vector<T> PuffCellVolume,
                           Array<vector<T>, 1 >& PuffCellConcentration);
    virtual void Forward();

    T ComputeInteractionCoefficient(Puff<T>* puff_alpha, Puff<T>* puff_beta);
    T ComputePuffOverlap(int alpha, int beta);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_GAUSSIANPUFFCHEMISTRY_HXX
#endif
