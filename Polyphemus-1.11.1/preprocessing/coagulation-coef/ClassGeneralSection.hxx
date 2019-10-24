
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_GENERAL_SECTION_HXX

namespace CoefficientRepartition
{
  class ClassParticle;

  class ClassGeneralSection : public ClassCoefficientRepartitionBase
  {
  public:
    // Nested vectors.
    //typedef CoefficientRepartition::double double;
    typedef CoefficientRepartition::Vector1I Vector1I;
    typedef CoefficientRepartition::Vector2I Vector2I;
    typedef CoefficientRepartition::Vector3I Vector3I;
    typedef CoefficientRepartition::Vector1T Vector1T;
    typedef CoefficientRepartition::Vector2T Vector2T;
    typedef CoefficientRepartition::Vector3T Vector3T;

  private:

    // Number of size sections.
    static int Nb_;

    // Number of composition sections.
    static int Nc_;

    // Bound section diameters.
    static Vector2T diameter_;

    // Bound section masses.
    static Vector2T mass_;

    // Vectors of composition index and fractions.
    static Vector3I fraction_index_;
    static Vector3T fraction_;

    // Index of size section.
    int size_bin_;

    // Index of composition bin.
    int composition_bin_;

  public:

    // Constructors.
    ClassGeneralSection();
    ClassGeneralSection(const int &size_bin, const int &fraction_bin);

    // Destructor.
    ~ClassGeneralSection();

    // Init.
    static void Init(const string type = "default");

    // Copy.
    void Copy(const ClassGeneralSection &gs);

    // Get methods.
    int GetNd() const;
    int GetSizeBin() const;
    int GetCompositionBin() const;
    double GetMassBoundary(const int &i) const;
    double GetDiameterBoundary(const int &i) const;
    void CollectFractionIndexBoundary(const int &i, Vector1I &fraction_index) const;
    void CollectFractionBoundary(const int &i, Vector1T &fraction) const;
    static int GetNdStatic(const int &c);
    static int GetNb();
    static int GetNc();
    static void CollectDiameter(Vector1T &diameter);
    static void CollectMass(Vector1T &mass);
    static void CollectFractionIndexBoundaryStatic(const int &c, const int &i, Vector1I &fraction_index);
    static void CollectFractionBoundaryStatic(const int &c, const int &i, Vector1T &fraction);

    // Randomly generate particles within section.
#ifndef SWIG
    inline void generate_random_particle(ClassParticle &p) const;
#endif
    void GenerateRandomParticle(ClassParticle &p) const;

    // Whether one particle is inside general section.
#ifndef SWIG
    inline bool has_particle(const ClassParticle &p) const;
#endif
    bool HasParticle(const ClassParticle &p) const;
  };
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_GENERAL_SECTION_HXX
#endif
