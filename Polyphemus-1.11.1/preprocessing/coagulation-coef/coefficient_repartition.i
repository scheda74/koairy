#ifdef WITH_MPI
%module(directors="1") coefficient_repartition_mpi_core
#else
%module(directors="1") coefficient_repartition_core
#endif
%{
#include "CoefficientRepartitionHeader.hxx"
  %}

%include "seldon.i"
%include <typemaps.i>
%include "std_iostream.i"
%include "std_string.i"
%include "std_vector.i"

namespace std
{
  %template(vector_int) vector<int>;
  %template(vector_str) vector<string>;
}

using namespace std;

%exception
{
  try
    {
      $action
        }
  catch(CoefficientRepartition::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Seldon::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(Ops::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "TalosHeader.hxx"
%include "SeldonHeader.hxx"

namespace CoefficientRepartition
{
  %rename(CoefficientRepartitionBase)  ClassCoefficientRepartitionBase;
  %rename(Particle)                    ClassParticle;
  %rename(GeneralSection)              ClassGeneralSection;
  %rename(CoefficientRepartition)      ClassCoefficientRepartition;

  %feature("director") ClassCoefficientRepartitionBase;
  %feature("director") ClassCoefficientRepartition;
  %feature("director") ClassParticle;
  %feature("director") ClassGeneralSection;
}

%include "CoefficientRepartitionHeader.hxx"
%include "ClassCoefficientRepartitionBase.hxx"
%include "ClassParticle.hxx"
%include "ClassGeneralSection.hxx"
%include "ClassCoefficientRepartition.hxx"
