
#ifndef COEFFICIENT_REPARTITION_FILE_ERROR_CXX

#include "Error.hxx"

namespace CoefficientRepartition
{
  Error::Error(string message = "") throw(): message_(message)
  {
  }


  Error::~Error() throw()
  {
  }


  string Error::What()
  {
    return "Error in CoefficientRepartition : " + message_;
  }


  void Error::CoutWhat()
  {
    cout << this->What() << endl;
  }
}

#define COEFFICIENT_REPARTITION_FILE_ERROR_CXX
#endif
