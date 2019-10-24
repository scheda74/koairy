//////////////
// INCLUDES //

#include <vector>
#include "species.h"
#include "solving.cxx"
using namespace soap;

void soap_main(double LWC, double RH, double Temperature, 
               double ionic, double chp, double& LWCorg,
               model_config* psoap_config,
               vector<species>* psurrogate, double& deltat,
               double DSD[], double csol[], double liquid[],
               int ns_aer, int neq, double q[], double qaero[], double qgas[]);

