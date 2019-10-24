#include <iostream>
using namespace std;

#define SELDONDATA_DEBUG_LEVEL_4
#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_ABORT
#define VERDANDI_WITH_ABORT
#define VERDANDI_DENSE
#include "Verdandi.hxx"
#include "method/OptimalInterpolation.cxx"
#include "observation/GroundNetworkObservationManager.cxx"

#include "Chimere.cxx"

//using namespace Polyphemus;


int main(int argc, char** argv)
{
  TRY;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  Verdandi::OptimalInterpolation < double, Polyphemus::Chimere,
                                   Polyphemus::GroundNetworkObservationManager<double> > driver(argv[1]);

  driver.Initialize();

  while (!driver.HasFinished())
    {
      driver.InitializeStep();
      driver.Forward();

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0)
        driver.Analyze();
    }

  END;

  return 0;
}
