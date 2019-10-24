
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_HXX

namespace CoefficientRepartition
{
  class ClassCoefficientRepartition : public ClassCoefficientRepartitionBase
  {
  public:
    // Nested vectors.
    //typedef CoefficientRepartition::double double;
    typedef CoefficientRepartition::Vector1I Vector1I;
    typedef CoefficientRepartition::Vector2I Vector2I;
    typedef CoefficientRepartition::Vector1T Vector1T;
    typedef CoefficientRepartition::Vector2T Vector2T;

  private:

    // Number of couples currently computed.
    int Ncompute_;

    // Monte Carlo number.
    int Nmc_;

    // Number of general sections.
    int Nsize_;

    // Vector of general sections.
    Vector<ClassGeneralSection, Vect_Full, NewAlloc<ClassGeneralSection> > general_section_;

    // Index of general sections.
    Vector2I index_first_, index_second_;

    // Repartition coefficient.
    Vector2T coefficient_;

  public:

    // Constructors.
    ClassCoefficientRepartition(const string type = "default",
                                const int Nmc = 100000);

    // Destructor.
    ~ClassCoefficientRepartition();

    // Get methods.
    int GetNcompute() const;
    int GetNsize() const;
    int GetNmc() const;
    ClassGeneralSection* GetGeneralSection(const int &i);
    void CollectIndexFirst(const int &i, Vector1I &index) const;
    void CollectIndexSecond(const int &i, Vector1I &index) const;
    void CollectCoefficient(const int &i, Vector1T &coefficient) const;

    // Clear.
    void Clear();

    // Compute repartition coefficients between general sections i1 and i2.
    void Compute(const int &i1, const int &i2);
    void ComputeAll();
#ifdef WITH_MPI
    void MPI_ComputeAll(const int recv_rank = 0);
#endif

    // Read and write coefficients.
#ifdef WITH_NETCDF
    void ReadNetCDF(const string &input_file);
    void WriteNetCDF(const string &output_file) const;
    void WriteTXT(const string &output_file) const;
    void WriteBIN(const string &output_file) const;
#endif
  };
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_COEFFICIENT_REPARTITION_HXX
#endif
