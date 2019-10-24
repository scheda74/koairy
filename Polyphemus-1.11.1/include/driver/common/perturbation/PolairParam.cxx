// Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
// Author(s): Mohamed Aissaoui, Vivien Mallet, Lin Wu
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


#ifndef POLYPHEMUS_FILE_PERTURBATION_POLAIRPARAM_CXX


#include "PolairParam.hxx"


namespace Polyphemus
{


  //! Default constructor.
  template <class T>
  PolairParam<T>::PolairParam():
    set_(false)
  {
  }


  //! Basic constructor.
  /*!
    This constructor just sets the parameter name.
    \param Name name of the parameter.
  */
  template <class T>
  PolairParam<T>::PolairParam(string name):
    name_(name), set_(false)
  {
  }


  //! Main constructor.
  /*!
    \param FieldName field name associated with the parameter.
    \param Name name of the parameter.
    \param Std standard deviation associated with the error on the parameter.
    \param Pdf type of the probability density function associated with the
    parameter.
    \param Param vector of pointers to parameters with which the current
    parameter is correlated.
    \param Correlation vector of correlations with other parameters.
    \param RefParam index of the reference parameter in 'Param'.
  */
  template <class T>
  PolairParam<T>::PolairParam(string FieldName, string Name, T Std,
                              string Pdf, vector<PolairParam<T>*> Param,
                              vector<T> Correlation, int RefParam):
    field_name_(FieldName), name_(Name), std_(Std), pdf_(Pdf), param_(Param),
    correlation_(Correlation), ref_param_(RefParam), set_(false)
  {
  }


  //! Destructor.
  template <class T>
  PolairParam<T>::~PolairParam()
  {
  }


  ////////////////////
  // ACCESS METHODS //
  ////////////////////


  //! Sets the field name associated with the parameter.
  /*!
    \param field_name field name associated with the parameter.
  */
  template <class T>
  void PolairParam<T>::SetFieldName(string field_name)
  {
    name_ = field_name;
  }


  //! Sets the parameter name.
  /*!
    \param Name name of the parameter.
  */
  template <class T>
  void PolairParam<T>::SetName(string name)
  {
    name_ = name;
  }


  //! Sets the type of the probability density function associated with the
  //! parameter.
  /*!
    \param pdf type of the probability density function associated with the
    parameter.
  */
  template <class T>
  void PolairParam<T>::SetPdf(string pdf)
  {
    pdf_ = pdf;
  }


  //! Sets the standard deviation associated with the error on the parameter.
  /*!
    \param std standard deviation associated with the error on the parameter.
  */
  template <class T>
  void PolairParam<T>::SetStd(T std)
  {
    std_ = std;
  }


  //! Sets the vector of pointers to parameters with which the current
  //! parameter is correlated.
  /*!
    \param param vector of pointers to parameters with which the current
    parameter is correlated.
  */
  template <class T>
  void PolairParam<T>::
  SetCorrelatedParam(const vector<PolairParam<T>*>& param)
  {
    param_ = param;
  }


  //! Sets the vector of correlations with other parameters.
  /*!
    \param correlation vector of correlations with other parameters.
  */
  template <class T>
  void PolairParam<T>::SetCorrelation(const vector<T>& correlation)
  {
    correlation_ = correlation;
  }


  //! Sets the index of the reference parameter in 'Param_'.
  /*!
    \param ref_param index of the reference parameter in 'Param_'.
  */
  template <class T>
  void PolairParam<T>::SetRefParam(int ref_param)
  {
    ref_param_ = ref_param;
  }


  //! Sets the random number which is used to modify the parameter field.
  /*!
    \param random_value the random number which is used to modify the field.
  */
  template <class T>
  void PolairParam<T>::SetRandomValue(T random_value)
  {
    random_value_ = random_value;
    set_ = true;
  }


  //! Returns the field name associated with the parameter
  /*!
    \return The field name associated with the parameter
  */
  template <class T>
  string PolairParam<T>::GetFieldName()
  {
    return field_name_;
  }


  //! Returns the parameter name.
  /*!
    \return The name of the parameter.
  */
  template <class T>
  string PolairParam<T>::GetName()
  {
    return name_;
  }


  //! Returns the type of the probability density function associated with the
  //! parameter.
  /*!
    \return The type of the probability density function associated with the
    parameter.
  */
  template <class T>
  string PolairParam<T>::GetPdf()
  {
    return pdf_;
  }


  //! Returns the standard deviation associated with the error on the
  //! parameter.
  /*!
    \return The standard deviation associated with the error on the parameter.
  */
  template <class T>
  T PolairParam<T>::GetStd()
  {
    return std_;
  }


  //! Returns the vector of pointers to parameters with which the current
  //! parameter is correlated.
  /*!
    \return The vector of pointers to parameters with which the current
    parameter is correlated.
  */
  template <class T>
  vector<PolairParam<T>*>& PolairParam<T>::GetCorrelatedParam()
  {
    return param_;
  }


  //! Returns the vector of correlations with other parameters.
  /*!
    \return The vector of correlations with other parameters.
  */
  template <class T>
  vector<T>& PolairParam<T>::GetCorrelation()
  {
    return correlation_;
  }


  //! Returns the index of the reference parameter in 'Param_'.
  /*!
    \return The index of the reference parameter in 'Param_'.
  */
  template <class T>
  int PolairParam<T>::GetRefParam()
  {
    return ref_param_;
  }


  //! Returns the random number which is used to modify the parameter field.
  /*!
    \return The the random number which is used to modify the parameter field.
  */
  template <class T>
  T PolairParam<T>::GetRandomValue()
  {
    return random_value_;
  }


  //! Is the random value set?
  /*!
    \return true is the random value was set, false otherwise.
  */
  template <class T>
  bool PolairParam<T>::IsSet()
  {
    return set_;
  }


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_POLAIRPARAM_CXX
#endif
