
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_PARTICLE_HXX

namespace CoefficientRepartition
{
  class ClassParticle : public ClassCoefficientRepartitionBase
  {
  public:
    // Nested vectors.
   // typedef CoefficientRepartition::double double;
    typedef CoefficientRepartition::Vector1T  Vector1T;

  private:

    // Particle mass.
    double mass_;

    // Vector of fraction species.
    Vector1T fraction_;

  public:

    // Allow direct access for these ClassGeneralSection method.
    friend void ClassGeneralSection::generate_random_particle(ClassParticle &p) const;
    friend bool ClassGeneralSection::has_particle(const ClassParticle &p) const;

    // Constructors.
    ClassParticle();
    ClassParticle(const double &mass, const Vector1T &fraction);

    // Destructor.
    ~ClassParticle();

    // Get methods.
    double GetDiameter() const;
    double GetMass() const;
    double GetFractionSpecies(const int &s) const;
    void CollectFraction(Vector1T &fraction) const;

    // Set methods.
    void SetData(const double &mass, const Vector1T &fraction);

    // Function coagulating two particles.
#ifndef SWIG
    inline ClassParticle coagulate(const ClassParticle &particle) const;
#endif
    ClassParticle Coagulate(const ClassParticle &particle) const;
  };
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_PARTICLE_HXX
#endif
