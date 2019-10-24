// Copyright (C) 2008-2009, INERIS - INRIA
// Author(s): Édouard Debry, Vivien Mallet
//
// This file is part of the air quality modeling system Polyphemus.
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


#ifndef POLYPHEMUS_FILE_MODELS_CHIMERE_HXX


#include "VerdandiHeader.hxx"
#include "BaseModel.cxx"
#include <mpi.h>


//////////////////////
// EXTERN FUNCTIONS //
//////////////////////


extern "C"
{
  void init_();
  void forward_();
  void initstep_();
  void cleanup_();
  void get_parameters_(int*, int*, int*, int*, int*, int*, int*, int*,
                       int*, int*, int*, int*, int*, int*, int*, int*);
  void get_species_(char*, int);
}

#include "chimere_pointers.hxx"

namespace Polyphemus
{


  /////////////
  // CHIMERE //
  /////////////


  /*! \brief 'Chimere' is an interface to the chemistry-transport model
    Chimere. */
  class Chimere: public Verdandi::VerdandiBase, public BaseModel<float>
  {
  public:
    typedef Verdandi::Vector<double> state;
    typedef Verdandi::Vector<double> state_error_variance_row;
    typedef Verdandi::Matrix<double> matrix_state_observation;

  protected:
    //! List of species as in Chimere.
    vector<string> species_list_chimere;

    //! Process rank in parallel computation.
    int rank;

    //! Number of parallel instances.
    int rank_size;

    //! Number of active species (gas and aerosols) as in Chimere.
    int Nspec;

    //! Total number of species.
    int Nspectot;

    //! Number of gas-chemistry reactions as in Chimere.
    int Nreac;

    //! Following integers values are defined in "src/chimere_params.f90.sed".
    //! Max number of rate constants.
    int Ntabmax;
    //! Max number of tabulated photolysis levels.
    int Nlevphotmax;
    //! Max number of photolysis reactions.
    int Nphotmax;
    //! Max number of tabulated zenith angle.
    int Ntabuzenmax;
    //! Number of prescribed species (O2, H2O...)
    int Nspresc;

    //! Number of components in aerosols as in Chimere.
    int Nkc;

    //! Number of bins as in Chimere.
    int Ms;

    //! Domain dimensions as in Chimere.
    int Nzonal, Nmerid, Nverti;

    //! Number of lateral boundary cells as in Chimere.
    int Nlatbound;

    //! List of aerosol bins.
    vector<string> bin_list;

    //! Bins bounds (in m).
    Array<double, 1> BinBound_aer;

    //! Number of species (gas and aerosol) with anthropogenic emissions.
    int Nemisa;

    //! List of species (gas and aerosol) with anthropogenic emissions.
    vector<string> species_list_emissions_anthropogenic;

    //! Number of species (gas and aerosol) with anthropogenic emissions.
    int Nemisb;

    //! List of species (gas and aerosol) with biogenic emissions.
    vector<string> species_list_emissions_biogenic;

    //! Number of level emissions as in Chimere.
    int Nlevemis;

    //! Number of Chimere families.
    int Nfam;

    //! List of families names as given in "chemprep/inputdata.XXXXXXXXX.X/FAMILIES"
    vector<string> species_list_families;

    //! Concentrations.
    Array<double, 4> Concentration_chimere;

    //! Maps of concentrations.
    map<string, Array<double, 3>* > Concentration_chimere_map;

    //! Maps of lateral boundaries.
    map<string, Array<float, 1>* > Lateral_boundary_chimere_map;

    //! Reaction rates.
    Array<double, 2> Reaction_rate_chimere;

    /*** Assimilation parameters ***/

    //! Balgovind length scale.
    double length_scale;
    //! State error variance.
    double state_error_variance;

  public:

    /*** Constructor and destructor ***/

    Chimere();
    Chimere(string config_file);
    void Construct(string config_file);
    virtual ~Chimere();

    /*** Initializations ***/

    virtual void Init();
    void Initialize(string config_file);
    virtual void InitStep();
    void InitializeStep();

    /*** Integration ***/

    virtual void Forward();

    /*** Finalize ***/

    bool HasFinished() const;
    void Finalize();

    /*** Configuration ***/

    virtual void ReadConfiguration();

    /*** Interface methods ***/

    void LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                            string name);
    void LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                            string name, int offset);
    void LinkFortranPointer(float* fortran_pointer, int N1, int N2, int N3,
                            string name);
    void LinkFortranPointer(float* fortran_pointer, int N1, int N2, int N3,
                            int N4, string name);
    void LinkFortranPointer(float* fortran_pointer, int N1, int N2, int N3,
                            string name, int offset);
    void LinkConcentration(double* fortran_pointer, int N1, int N2, int N3,
                           int N4);
    void LinkFortranPointer(float* fortran_pointer, int N1, int N2, int N3,
                            int N4, int N5, string name);

    /*** Access methods ***/

    double GetTime();
    Date GetDate(double time);
    void GetState(Array<float, 1>& state_vector);
    void GetState(state& state_vector);
    void SetState(Array<float, 1> state_vector);
    void SetState(const state& state_vector);
    int GetNstate() const;
    Array<double, 3>& GetConcentrationChimere(string field);
    Array<float, 1>& GetLateralBoundaryChimere(string field);

    void GetStateErrorVarianceRow(int row, state_error_variance_row& P_row);

    string GetName() const;
    void Message(string message);
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CHIMERE_HXX
#endif
