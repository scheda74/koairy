
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_BASE_HXX

namespace CoefficientRepartition
{
  class ClassCoefficientRepartitionBase
  {
  public:
    //typedef CoefficientRepartition::double double;

  protected:

    // MPI variables.
    static int rank_;
    static int Nrank_;

    // Configuration file.
    static string configuration_file_;

    // Number of species.
    static int Ns_;

    // Particle density.
    static double density_;

  public:
    // Initialize.
    static void Init(const string &configuration_file, string config_type = "default");

#ifdef WITH_MPI
    static void MPI_Init();
#endif

    // Constructors.
    ClassCoefficientRepartitionBase();

    // Destructor.
    ~ClassCoefficientRepartitionBase();

    // Get methods.
    static int GetRank();
    static int GetNrank();
    static int GetNs();
    static double GetDensity();
  };
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_BASE_HXX
#endif
