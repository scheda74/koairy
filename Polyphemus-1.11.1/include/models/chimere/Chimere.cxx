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


#ifndef POLYPHEMUS_FILE_MODELS_CHIMERE_CXX


#include "Verdandi.hxx"

#include "Chimere.hxx"

#define WORKTAG 2

namespace Polyphemus
{


  Chimere::Chimere()
  {
  }


  Chimere::Chimere(string config_file)
  {
    Construct(config_file);
  }


  void Chimere::Construct(string config_file)
  {
    init_();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);

    BaseModel<float>::Construct(config_file);

    Nspec = 0;
    Nspectot = 0;
    Nreac = 0;
    Nkc = 0;
    Ms = 0;
    Nzonal = 0;
    Nmerid = 0;
    Nverti = 0;
    Nlatbound = 0;
    Nemisa = 0;
    Nemisb = 0;
    Nlevemis = 0;
    Nfam = 0;
    Nspresc = 4;
    Ntabmax = 22;
    Nlevphotmax = 50;
    Nphotmax = 50;
    Ntabuzenmax = 20;

    // Assimilation-related parameters.
    config.SetSection("[data_assimilation]");
    config.PeekValue("Length_scale", length_scale);
    config.PeekValue("State_error_variance", state_error_variance);
  }


  Chimere::~Chimere()
  {
    Finalize();
  }


  void Chimere::Init()
  {
    BaseModel<float>::Init();

    Nspec = this->Ns;
    if (Nkc * Ms != 1)
      Nspec += Nkc * Ms;
    Nspectot = Nspec + Nspresc + Nfam;
    Nzonal = this->Nx;
    Nmerid = this->Ny;
    Nverti = this->Nz;
    Nlatbound = 2 * (Nmerid + Nzonal) * Nverti;

    // Get and check parameters.
    int Nspec_, Nreac_, Nkc_, Ms_, Nzonal_, Nmerid_,
      Nverti_, Nemisa_, Nemisb_, Nlevemis_, Nfam_;

    get_parameters_(&Nspresc, &Ntabmax, &Nlevphotmax, &Nphotmax, &Ntabuzenmax,
                    &Nspec_, &Nreac_, &Nkc_, &Ms_, &Nzonal_, &Nmerid_,
                    &Nverti_, &Nemisa_, &Nemisb_, &Nlevemis_, &Nfam_);

    if (Nspec != Nspec_)
      throw string("Nspec (" + to_str(Nspec) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nspec_) + ")");
    if (Nreac != Nreac_)
      throw string("Nreac (" + to_str(Nreac) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nreac_) + ")");
    if (Nkc != Nkc_)
      throw string("Nkc (" + to_str(Nkc) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nkc_) + ")");
    if (Ms != Ms_)
      throw string("Ms (" + to_str(Ms) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Ms_) + ")");
    if (Nzonal != Nzonal_)
      throw string("Nzonal (" + to_str(Nzonal) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nzonal_) + ")");
    if (Nmerid != Nmerid_)
      throw string("Nmerid (" + to_str(Nmerid) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nmerid_) + ")");
    if (Nverti != Nverti_)
      throw string("Nverti (" + to_str(Nverti) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nverti_) + ")");
    if (Nemisa != Nemisa_)
      throw string("Nemisa (" + to_str(Nemisa) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nemisa_) + ")");
    if (Nemisb != Nemisb_)
      throw string("Nemisb (" + to_str(Nemisb) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nemisb_) + ")");
    if (Nlevemis != Nlevemis_)
      throw string("Nlevemis (" + to_str(Nlevemis) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nlevemis_) + ")");
    if (Nfam != Nfam_)
      throw string("Nfam (" + to_str(Nfam) + ") from configuration does not match "
                   + "Chimere running value (" + to_str(Nfam_) + ")");

    if (rank == 0)
      {
        int len = 3000;
        char species_list_tmp[len];
        for (int pos = 0; pos < len; pos++)
          species_list_tmp[pos] = ' ';
        get_species_(species_list_tmp, len);
        string str_tmp(species_list_tmp);
        species_list_chimere = split(str_tmp, " ");
      }

    Range all = Range::all();
    string name;

    if (rank == 0)
      {
        LinkConcentration(chimere_common_MP_conc,
                          Nspectot, Nzonal + 6, Nmerid + 6, Nverti + 1);
        for (int i = 0; i < Nspectot; i++)
          {
            name = species_list_chimere[i];
            Concentration_chimere_map[name] = new Array<double, 3>;
            Concentration_chimere_map[name]
              ->reference(Concentration_chimere(all, all, all, i));
          }
        LinkFortranPointer(chimere_common_MP_airm,
                           Nzonal, Nmerid, Nverti,
                           "AirDensity", 1);
        LinkFortranPointer(chimere_common_MP_winz,
                           Nzonal, Nmerid, Nverti,
                           "ZonalWind", 1);
        LinkFortranPointer(chimere_common_MP_winm,
                           Nzonal, Nmerid, Nverti,
                           "MeridionalWind", 1);
        LinkFortranPointer(chimere_common_MP_hlay,
                           Nzonal, Nmerid, Nverti,
                           "LayerThickness", 1);
        LinkFortranPointer(chimere_common_MP_sphu,
                           Nzonal, Nmerid, Nverti,
                           "SpecificHumidity", 1);
        LinkFortranPointer(chimere_common_MP_temp,
                           Nzonal, Nmerid, Nverti,
                           "Temperature", 1);
        LinkFortranPointer(chimere_common_MP_kzzz,
                           Nzonal, Nmerid, Nverti,
                           "TopLayerKz", 1);
        LinkFortranPointer(chimere_common_MP_clwc,
                           Nzonal, Nmerid, Nverti,
                           "LiquidWaterContent", 1);
        LinkFortranPointer(chimere_common_MP_hght,
                           Nzonal, Nmerid,
                           "MixingHeight", 1);
        LinkFortranPointer(chimere_common_MP_atte,
                           Nzonal, Nmerid,
                           "Attenuation", 1);
        LinkFortranPointer(chimere_common_MP_tem2,
                           Nzonal, Nmerid,
                           "Temperature2m", 1);
        LinkFortranPointer(chimere_common_MP_usta,
                           Nzonal, Nmerid,
                           "FrictionVelocity", 1);
        LinkFortranPointer(chimere_common_MP_aerr,
                           Nzonal, Nmerid,
                           "AerodynamicResistance", 1);
        LinkFortranPointer(chimere_common_MP_obuk,
                           Nzonal, Nmerid,
                           "LMO", 1);
        LinkFortranPointer(chimere_common_MP_wsta,
                           Nzonal, Nmerid,
                           "ConvectiveVelocity", 1);
        LinkFortranPointer(chimere_common_MP_topc,
                           Nzonal, Nmerid,
                           "TotalPrecipitation", 1);
        LinkFortranPointer(chimere_common_MP_sreh,
                           Nzonal, Nmerid,
                           "SurfaceRelativeHumidity", 1);
        LinkFortranPointer(chimere_common_MP_emisaloc,
                           Nemisa, Nzonal, Nmerid, Nlevemis,
                           "AnthropogenicEmission");
        for (int i = 0; i < Nemisa; i++)
          {
            name = "AnthropogenicEmission_" + species_list_emissions_anthropogenic[i];
            this->A3_map[name] = new Array<float, 3>;
            this->A3_map[name]
              ->reference(this->A4("AnthropogenicEmission")(all, all, all, i));
          }
        LinkFortranPointer(chimere_common_MP_emisb,
                           Nemisb, Nzonal, Nmerid,
                           "BiogenicEmission", 1);
        for (int i = 0; i < Nemisb; i++)
          {
            name = "BiogenicEmission_" + species_list_emissions_biogenic[i];
            this->A2_map[name] = new Array<float, 2>;
            this->A2_map[name]
              ->reference(this->A3("BiogenicEmission")(all, all, i));
          }
        Reaction_rate_chimere
          .reference(Array<double, 2>(chimere_common_MP_tabrate + Ntabmax + 1,
                                      shape(Nreac, Ntabmax),  neverDeleteData));
        LinkFortranPointer(chimere_common_MP_photoj,
                           Ntabuzenmax, Nlevphotmax, Nphotmax,
                           "PhotolysisRate");
        LinkFortranPointer(chimere_common_MP_boundlat,
                           Nlatbound, Nspec,
                           "LateralBoundary", 1);
        for (int i = 0; i < Nspec; i++)
          {
            name = species_list_chimere[i];
            Lateral_boundary_chimere_map[name] = new Array<float, 1>;
            Lateral_boundary_chimere_map[name]
              ->reference(this->A2("LateralBoundary")(i, all));
          }
        LinkFortranPointer(chimere_common_MP_boundtop,
                           Nzonal, Nmerid, Nspec,
                           "TopBoundary", 1);
        for (int i = 0; i < Nspec; i++)
          {
            name = "TopBoundary_" + species_list_chimere[i];
            this->A2_map[name] = new Array<float, 2>;
            this->A2_map[name]
              ->reference(this->A3("TopBoundary")(i, all, all));
          }
      }

    // Initialization both in master and slaves.
    LinkFortranPointer(chimere_perturbation_MP_perturbation_depoloc,
                       Nspec, 1, "PerturbationDepositionVelocity");
    A2("PerturbationDepositionVelocity") = 1.0;

    LinkFortranPointer(chimere_perturbation_MP_perturbation_rate,
                       Nreac, 1, 1, "PerturbationReactionRate");
    A3("PerturbationReactionRate") = 1.0;
    for (int i = 0; i < Nreac; i++)
      {
        name = "PerturbationReactionRate_" + to_str(i + 1);
        this->A2_map[name] = new Array<float, 2>;
        this->A2_map[name]
          ->reference(this->A3("PerturbationReactionRate")(all, all, i));
      }
  }


  void Chimere::Initialize(string config_file)
  {
    Construct(config_file);
    Init();
  }


  void Chimere::InitStep()
  {
    BaseModel<float>::InitStep();
    initstep_();
    A2("PerturbationDepositionVelocity") = 1.0;
    A3("PerturbationReactionRate") = 1.0;
  }


  void Chimere::InitializeStep()
  {
    InitStep();
  }


  void Chimere::Forward()
  {
    int buf_size = 1000;
    float* buf = new float[buf_size];
    for (int i = 0; i < buf_size; i++)
      buf[i] = 0.0;

    if (rank == 0)
      {
        for (int i = 0; i < Nspec; i++)
          buf[i] = A2("PerturbationDepositionVelocity")(0, i);
        for (int i = 0; i < Nreac; i++)
          buf[Nspec + i] = A3("PerturbationReactionRate")(0, 0, i);

        for (int i = 1; i < rank_size; i++)
          MPI_Send(buf, buf_size, MPI_FLOAT, i, WORKTAG, MPI_COMM_WORLD);
      }
    else
      {
        MPI_Status status;
        MPI_Recv(buf, buf_size, MPI_FLOAT, 0, WORKTAG, MPI_COMM_WORLD, &status);
        for (int i = 0; i < Nspec; i++)
          A2("PerturbationDepositionVelocity")(0, i) = buf[i];
        for (int i = 0; i < Nreac; i++)
          A3("PerturbationReactionRate")(0, 0, i) = buf[Nspec + i];
      }

    BaseModel<float>::Forward();
    forward_();
    this->AddTime(this->Delta_t);
    this->step++;
  }


  bool Chimere::HasFinished() const
  {
    return this->step >= this->Nt;
  }


  void Chimere::Finalize()
  {
    cleanup_();
  }


  void Chimere::ReadConfiguration()
  {
    BaseModel<float>::ReadConfiguration();

    /*** Species aerosol ***/

    // Opens the file that describes species.
    ConfigStream species_stream(this->file_species);
    // Section "[species]" contains all species names.
    species_stream.SetSection("[aerosol_species]");
    while (!species_stream.IsEmpty())
      this->species_list_aer.push_back(species_stream.GetElement());
    this->Ns_aer = int(this->species_list_aer.size());
    if (this->Ns_aer == 0)
      this->Ns_aer = 1;
    Nkc = this->Ns_aer;

    // Reads bin bounds.
    if (Nkc != 1)
      {
        this->config.SetSection("[domain]");
        this->config.Find("Bin_bounds");
        bin_list = split(this->config.GetLine());
        this->Nbin_aer = int(bin_list.size()) - 1;

        BinBound_aer.resize(this->Nbin_aer + 1);
        // Reads bin bounds in micrometers and converts it to meters.
        for (int i = 0; i < this->Nbin_aer + 1; i++)
          BinBound_aer(i) = 1.e-6 * convert<double>(bin_list[i]);
      }
    else
      this->Nbin_aer = 1;
    Ms = this->Nbin_aer;

    BinBound_aer.resize(this->Nbin_aer + 1);
    // Reads the bin bounds in micrometers and converts it to meters.
    for (int i = 0; i < this->Nbin_aer + 1; i++)
      BinBound_aer(i) = 1.e-6 * convert<double>(bin_list[i]);

    // Reads the number of reactions.
    this->config.PeekValue("Nreac", Nreac);

    // Reads the path to the data configuration file.
    this->config.SetSection("[data]");
    string data_description_file = this->config.PeekValue("Data_description");

    // Opens the configuration file for input data.
    ConfigStream data_description_stream(data_description_file);

    // Emissions.
    data_description_stream.SetSection("[emissions]");
    data_description_stream.Find("anthropogenic");
    species_list_emissions_anthropogenic = split(data_description_stream.GetLine());
    Nemisa = int(species_list_emissions_anthropogenic.size());

    data_description_stream.Find("biogenic");
    species_list_emissions_biogenic = split(data_description_stream.GetLine());
    Nemisb = int(species_list_emissions_biogenic.size());

    data_description_stream.PeekValue("Nz", Nlevemis);

    // Families.
    data_description_stream.SetSection("[families]");
    while (!data_description_stream.IsEmpty())
      species_list_families.push_back(data_description_stream.GetElement());
    Nfam = int(species_list_families.size());
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   string name)
  {
    this->A2_map[name] = new Array<float, 2>;
    this->A2_map[name]
      ->reference(Array<float, 2>(fortran_pointer + N1 + 1,
                                  shape(N2, N1), neverDeleteData));
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   string name, int offset)
  {
    this->A2_map[name] = new Array<float, 2>;
    this->A2_map[name]
      ->reference(Array<float, 2>(fortran_pointer + N2 * N1 * (offset + 1) + N1 + 1,
                                  shape(N2, N1), neverDeleteData));
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   int N3, string name)
  {
    this->A3_map[name] = new Array<float, 3>;
    this->A3_map[name]
      ->reference(Array<float, 3>(fortran_pointer + N2 * N1 + N1 + 1,
                                  shape(N3, N2, N1), neverDeleteData));
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   int N3, string name, int offset)
  {
    this->A3_map[name] = new Array<float, 3>;
    this->A3_map[name]
      ->reference(Array<float, 3>(fortran_pointer + N3 * N2 * N1 * (offset + 1)
                                  + N2 * N1 + N1 + 1,
                                  shape(N3, N2, N1), neverDeleteData));
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   int N3, int N4, string name)
  {
    this->A4_map[name] = new Array<float, 4>;
    this->A4_map[name]
      ->reference(Array<float, 4>(fortran_pointer + N3 * N2 * N1
                                  + N2 * N1 + N1 + 1,
                                  shape(N4, N3, N2, N1), neverDeleteData));
  }


  void Chimere::LinkConcentration(double* fortran_pointer, int N1, int N2,
                                  int N3, int N4)
  {
    Concentration_chimere
      .reference(Array<double, 4>(fortran_pointer + N3 * N2 * N1
                                  + N2 * N1 + N1 + 1,
                                  shape(N4, N3, N2, N1), neverDeleteData));
  }


  void Chimere::LinkFortranPointer(float* fortran_pointer, int N1, int N2,
                                   int N3, int N4, int N5, string name)
  {
    this->A5_map[name] = new Array<float, 5>;
    this->A5_map[name]
      ->reference(Array<float, 5>(fortran_pointer + N4 * N3 * N2 * N1
                                  + N3 * N2 * N1 + N2 * N1 + N1 + 1,
                                  shape(N5, N4, N3, N2, N1),
                                  neverDeleteData));
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  double Chimere::GetTime()
  {
    return this->GetCurrentTime();
  }


  Date Chimere::GetDate(double time)
  {
    Date out;
    out = this->Date_min;
    out.AddSeconds(time);
    return out;
  }


  void Chimere::GetState(Array<float, 1>& state_vector)
  {
    state_vector = Array<float, 1>();
  }


  void Chimere::GetState(state& state_vector)
  {
    // Avogadro number S.I.
    double avogadro = 6.022045e23;

    // Ozone molar mass in µg/mol.
    double molar_mass = 48.0e6;

    // Conversion factor from molecules/cm^3 to µg/m^3.
    double conversion_factor = 1.e6 * molar_mass / avogadro;

    state_vector.Reallocate(this->Nx * this->Ny);
    int position = 0;
    for (int j = 0; j < this->Ny; j++)
      for (int i = 0; i < this->Nx; i++)
        state_vector(position++) = GetConcentrationChimere("O3")(0, j, i);

    state_vector *= conversion_factor;
  }


  void Chimere::SetState(Array<float, 1> state_vector)
  {
  }


  void Chimere::SetState(const state& state_vector)
  {
    // Avogadro number S.I.
    double avogadro = 6.022045e23;

    // Ozone molar mass in µg/mol.
    double molar_mass = 48.0e6;

    // Conversion factor from µg/m^3 to molecules/cm^3.
    double conversion_factor = 1.e-6 / molar_mass * avogadro;

    int position = 0;
    for (int j = 0; j < this->Ny; j++)
      for (int i = 0; i < this->Nx; i++)
        GetConcentrationChimere("O3")(0, j, i) = state_vector(position++)
          * conversion_factor;
  }


  int Chimere::GetNstate() const
  {
    return this->Ny * this->Nx;
  }


  Array<double, 3>& Chimere::GetConcentrationChimere(string field)
  {
    if (Concentration_chimere_map.find(field) == Concentration_chimere_map.end())
      throw string("Field \"") + field + string("\" cannot be found in Chimere.");
    else
      return *Concentration_chimere_map[field];
  }


  Array<float, 1>& Chimere::GetLateralBoundaryChimere(string field)
  {
    if (Lateral_boundary_chimere_map.find(field) == Lateral_boundary_chimere_map.end())
      throw string("Field \"") + field + string("\" cannot be found in Chimere.");
    else
      return *Lateral_boundary_chimere_map[field];
  }


  //! Computes a row of the variance of the state error.
  /*!
    \param[in] row row index.
    \param[out] P_row the row with index \a row in the state error variance.
  */
  void Chimere
  ::GetStateErrorVarianceRow(int row, state_error_variance_row& P_row)
  {
    Seldon::Vector<int> position;
    Seldon::Vector<int> N(2);
    N(0) = Ny;
    N(1) = Nx;
    Verdandi::get_position(row, N, position);
    int row_y = position(0);
    int row_x = position(1);

    P_row.Reallocate(Nx * Ny);
    for (int j = 0; j < Ny; j++)
      for (int i = 0; i < Nx; i++)
        {
          double ratio, distance(0), tmp;
          tmp = double(row_x - i) * Delta_x;
          distance += tmp * tmp;
          tmp = double(row_y - j) * Delta_y;
          distance += tmp * tmp;
          ratio = sqrt(distance) / length_scale;
          P_row(i + Nx * j) = state_error_variance * (1. + ratio)
            * exp(-ratio);
        }
  }


  //! Returns the name of the class, that is, "Chimere".
  /*!
    \return The name of the class, that is, "Chimere".
  */
  string Chimere::GetName() const
  {
    return "Chimere";
  }


  //! Receives and handles a message.
  /*
    \param[in] message the received message.
  */
  void Chimere::Message(string message)
  {
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_CHIMERE_CXX
#endif
