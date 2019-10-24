#ifndef COEFFICIENT_REPARTITION_FILE_ERROR_HXX

namespace CoefficientRepartition
{
  class Error: public std::exception
  {
  private:
    string message_;
  public:
    Error(string message) throw();

    // Destructor.
    virtual ~Error() throw();
    
    virtual string What();
    void CoutWhat();
  };
}

#define COEFFICIENT_REPARTITION_FILE_ERROR_HXX
#endif
