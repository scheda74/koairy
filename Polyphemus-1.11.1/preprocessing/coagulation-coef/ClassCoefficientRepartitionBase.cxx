
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_BASE_CXX

namespace CoefficientRepartition
{
  // Initialize.
  void ClassCoefficientRepartitionBase::Init(const string &configuration_file, string config_type)
  {
    configuration_file_ = configuration_file;

    Ops::Ops ops(configuration_file_);
    Ns_ = ops.Get<int>("Number_species", "", Ns_);

    // Density converted from g/cm3  to µg/µm3.
    density_ = ops.Get<double>("Particle_density", "", density_) * 1.e-6;
    ops.Close();

#ifdef WITH_MPI
    rank_ = MPI::COMM_WORLD.Get_rank();
    Nrank_ = MPI::COMM_WORLD.Get_size();
#else
    rank_ = 0;
    Nrank_ = 1;
#endif

    ClassGeneralSection::Init(config_type);
  }


#ifdef WITH_MPI
  void ClassCoefficientRepartitionBase::MPI_Init()
  {
    MPI::Init();
  }
#endif

  // Constructors.
  ClassCoefficientRepartitionBase::ClassCoefficientRepartitionBase()
  {
    return;
  }


  // Destructor.
  ClassCoefficientRepartitionBase::~ClassCoefficientRepartitionBase()
  {
    return;
  }


  // Get methods.
  int ClassCoefficientRepartitionBase::GetRank()
  {
    return rank_;
  }


  int ClassCoefficientRepartitionBase::GetNrank()
  {
    return Nrank_;
  }


  int ClassCoefficientRepartitionBase::GetNs()
  {
    return Ns_;
  }


  double ClassCoefficientRepartitionBase::GetDensity()
  {
    return density_;
  }
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_BASE_CXX
#endif
