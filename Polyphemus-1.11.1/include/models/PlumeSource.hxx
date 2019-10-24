// Copyright (C) 2006-2012, ENPC - INRIA - EDF R&D
// Author(s): Hadjira Foudhil, Vivien Mallet, Régis Briant
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

// This file is part of a Gaussian plume model for Polyphemus.


#ifndef POLYPHEMUS_FILE_MODELS_PLUMESOURCE_HXX


#include <list>


namespace Polyphemus
{


  using namespace std;


  //////////////////
  // SOURCE PLUME //
  //////////////////


  //! This class stores a source description.
  template<class T>
  class PlumeSource
  {
  protected:

    //! Rate.
    T rate_;
    //! Efflux speed of gases.
    T velocity_;
    //! Temperature.
    T temperature_;
    //! Diameter.
    T diameter_;
    //! Abscissa.
    T x_;
    //! Ordinate.
    T y_;
    //! Height.
    T z_;
    T effective_height_;
    //! Height of the fraction of plume above BL.
    T effective_height_above_;
    //! Fraction of the plume above BL.
    T penetration_factor_;
    //! Square of the initial horizontal plume spread.
    T sigma_y_2_;
    //! Square of the initial vertical plume spread.
    T sigma_z_2_;
    //! Species index.
    int species_index_;
    //! Type
    string type_;
    //! Identifier of the source.
    int id_source_;


  public:

    /*** Constructor ***/

    PlumeSource(T rate, T velocity, T temperature, T diameter,
                T x, T y, T z, int species_index, int id_source = 0);

    /*** Methods ***/

    virtual ~PlumeSource();
    T GetRate() const;
    void SetRate(T rate);
    void SetPenetrationFactor(T P);
    T GetPenetrationFactor() const;
    T GetVelocity() const;
    void SetVelocity(T velocity);
    T GetTemperature() const;
    void SetTemperature(T temperature);
    T GetDiameter() const;
    void SetDiameter(T diameter);
    T GetX() const;
    void SetX(T x);
    T GetY() const;
    void SetY(T y);
    T GetZ() const;
    void SetZ(T z);
    T GetHeight() const;
    T GetHeightAboveBL() const;
    void SetHeight(T effective_height);
    void SetHeightAboveBL(T effective_height_above);
    T GetSigma_y() const;
    T GetSigma_z() const;
    void SetSigma_y_2(T sigma_y);
    void SetSigma_z_2(T sigma_z);
    T GetSigma_y_2() const;
    T GetSigma_z_2() const;
    int GetSpeciesIndex() const;
    virtual T GetX2() const;
    virtual T GetY2() const;
    virtual T GetZ2() const;
    virtual void SetX2(T X_2);
    virtual void SetY2(T Y_2);
    virtual void SetZ2(T Z_2);
    string GetEmissionType() const;
    virtual T GetCombinationFactor() const;
    virtual void SetCombinationFactor(T combination_factor);
    virtual T GetWidth() const;
    virtual int GetIdSource() const;
    virtual int GetIdSection() const;
    virtual T GetVehicleVelocity() const;
    virtual T GetArea() const;
    virtual T GetDensity() const;
  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_PLUMESOURCE_HXX
#endif
