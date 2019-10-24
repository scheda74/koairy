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


#ifndef POLYPHEMUS_FILE_PERTURBATION_POLAIRPARAM_HXX


namespace Polyphemus
{


  /////////////////
  // POLAIRPARAM //
  /////////////////


  //! 'PolairParam' stores various data associated with a parameter.
  template<class T>
  class PolairParam
  {

  protected:

    //! Field name associated with the parameter.
    string field_name_;
    //! Parameter name.
    string name_;
    //! Standard deviation associated with the error of the parameter.
    T std_;
    /*! \brief Type of the probability density function associated with the
      parameter.
    */
    string pdf_;
    /*! \brief Vector of pointers to parameters with which the current
      parameter is correlated.
    */
    vector<PolairParam<T>*> param_;
    //! Vector of correlations with other parameters.
    vector<T> correlation_;
    //! Index of the reference parameter in 'param_'.
    int ref_param_;

    //! Random number which is used to modify the field.
    T random_value_;
    //! Has the random value been set?
    bool set_;

  public:

    /*** Constructor ***/

    PolairParam();
    PolairParam(string name);
    PolairParam(string FieldName, string Name, T Std, string Pdf,
                vector<PolairParam<T>*> Param = vector<PolairParam<T>*>(),
                vector<T> Correlation = vector<T>(), int RefParam = -1);

    virtual ~PolairParam();

    /***  Access methods ***/

    void SetFieldName(string field_name);
    void SetName(string name);
    void SetPdf(string pdf);
    void SetStd(T std);
    void SetCorrelatedParam(const vector<PolairParam<T>*>& param);
    void SetCorrelation(const vector<T>& Correlation);
    void SetRefParam(int RefParam);
    void SetRandomValue(T random_value);

    string GetFieldName();
    string GetName();
    string GetPdf();
    T GetStd();
    vector<PolairParam<T>*>& GetCorrelatedParam();
    vector<T>& GetCorrelation();
    int GetRefParam();
    T GetRandomValue();

    bool IsSet();

  };


} // namespace Polyphemus.


#define POLYPHEMUS_FILE_PERTURBATION_POLAIRPARAM_HXX
#endif
