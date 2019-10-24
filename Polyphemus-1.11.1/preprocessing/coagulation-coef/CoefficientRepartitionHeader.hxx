
#ifndef COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_HEADER_HXX

namespace std
{
}

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <cstdlib>
#include <string>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <map>
#include <vector>
#include <utility>

#include "talos/TalosHeader.hxx"
using namespace Talos;

#include "seldon/SeldonHeader.hxx"
using namespace Seldon;

#define OPS_WITH_EXCEPTION
#include "ops/OpsHeader.hxx"

// MPI.
#ifdef WITH_MPI
#define MPI2 1
#define COEFFICIENT_REPARTITION_MPI_REAL MPI::DOUBLE
#define MPI_INCLUDED
#include <mpi.h>
#else
#define MPI2 0
#endif

// NetCDF library.
#ifdef WITH_NETCDF
#define NETCDF 1
#include <netcdfcpp.h>
#define COEFFICIENT_REPARTITION_NETCDF_REAL ncDouble
#else
#define NETCDF 0
#endif

namespace CoefficientRepartition
{
  using namespace std;

#define PI_6     0.52359877559829882
#define INV_PI_6 1.909859317102744
#define FRAC3    0.33333333333333333

  typedef Vector<int> Vector1I;
  typedef Vector<Vector1I, Vect_Full, NewAlloc<Vector1I > > Vector2I;
  typedef Vector<Vector2I, Vect_Full, NewAlloc<Vector2I > > Vector3I;
  typedef Vector<double> Vector1T;
  typedef Vector<Vector1T, Vect_Full, NewAlloc<Vector1T > > Vector2T;
  typedef Vector<Vector2T, Vect_Full, NewAlloc<Vector2T > > Vector3T;

#ifndef SWIG
  using Talos::to_str;
  using Talos::to_num;
#endif
}

#include "Error.hxx"
#include "ClassCoefficientRepartitionBase.hxx"
#include "ClassGeneralSection.hxx"
#include "ClassParticle.hxx"
#include "ClassCoefficientRepartition.hxx"

#ifdef TRY
#undef TRY
#endif
#define TRY try                                         \
    {

#ifdef END
#undef END
#endif
#define END                                             \
  }                                                     \
    catch (CoefficientRepartition::Error& err)          \
      {                                                 \
        err.CoutWhat();                                 \
        return 1;                                       \
      }                                                 \
    catch (Ops::Error& err)	                        \
      {                                                 \
        err.CoutWhat();                                 \
        return 1;                                       \
      }                                                 \
    catch (Seldon::Error& err)                          \
      {                                                 \
        err.CoutWhat();                                 \
        return 1;                                       \
      }                                                 \
    catch (std::exception& err)                         \
      {                                                 \
        cout << "C++ exception: "                       \
             << err.what() << endl;                     \
        return 1;                                       \
      }                                                 \
    catch (std::string& str)                            \
      {                                                 \
        cout << str << endl;                            \
        return 1;                                       \
      }                                                 \
    catch (const char* str)                             \
      {                                                 \
        cout << str << endl;                            \
        return 1;                                       \
      }                                                 \
    catch(...)                                          \
      {                                                 \
        cout << "Unknown exception..." <<endl;          \
        return 1;                                       \
      }

#define COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_HEADER_HXX
#endif
