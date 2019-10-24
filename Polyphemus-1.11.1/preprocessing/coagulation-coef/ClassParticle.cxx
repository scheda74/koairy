
#ifndef COEFFICIENT_REPARTITION_FILE_CLASS_PARTICLE_CXX

#include "ClassParticle.hxx"

namespace CoefficientRepartition
{
  // Constructors.
  ClassParticle::ClassParticle() : ClassCoefficientRepartitionBase()
  {
    // Default is a one micrometer uniform composition particle.
    mass_ = this->density_ * PI_6;

    fraction_.Reallocate(this->Ns_);
    fraction_.Fill(double(1) / double(this->Ns_));

    return;
  }


  ClassParticle::ClassParticle(const double &mass, const Vector1T &fraction)
    : ClassCoefficientRepartitionBase()
  {
    mass_ = mass;

    if (fraction.GetSize() != this->Ns_)
      throw CoefficientRepartition::Error("Fraction vector size is not equal to the number of species.");

    if (Norm1(fraction) != double(1))
      throw CoefficientRepartition::Error("Fraction sum is not equal to unity.");

    fraction_ = fraction;

    return;
  }


  // Destructor.
  ClassParticle::~ClassParticle()
  {
    return;
  }


  // Get methods.
  double ClassParticle::GetDiameter() const
  {
    return pow(mass_ / this->density_ * INV_PI_6, FRAC3);
  }


  double ClassParticle::GetMass() const
  {
    return mass_;
  }


  double ClassParticle::GetFractionSpecies(const int &s) const
  {
    return fraction_(s);
  }


  void ClassParticle::CollectFraction(Vector1T &fraction) const
  {
    fraction = fraction_;
  }


  // Set methods.
  void ClassParticle::SetData(const double &mass, const Vector1T &fraction)
  {
    mass_ = mass;
    fraction_ = fraction;
  }


  // Function coagulating two particles.
  ClassParticle ClassParticle::Coagulate(const ClassParticle &particle) const
  {
    return coagulate(particle);
  }


  ClassParticle ClassParticle::coagulate(const ClassParticle &particle) const
  {
    ClassParticle coag;

    coag.mass_ = mass_ + particle.mass_;
    //cout<<"mass 1="<<mass_<<"mass 2="<<particle.mass_<<endl;
    double tmp_fraction(0);
    for (int s = 0; s < this->Ns_; s++)
    {
      coag.fraction_(s) = (fraction_(s) * mass_ + particle.fraction_(s) *  particle.mass_) / coag.mass_;
      tmp_fraction+=coag.fraction_(s);
      //cout<<"coag.fraction_("<<s<<")="<<coag.fraction_(s)<<endl;
      if(coag.fraction_(s)<0.0)
      {
	cout<<"error, coagulation: "<<" frac1="<<fraction_(s)<<" frac2="<<particle.fraction_(s)<<" "<<coag.mass_<<endl;
	int i;
	cin>>i; 	
      }
    }
    //if(tmp_fraction!=1.0)
    //{
	//  cout<<"error, coagulation "<<tmp_fraction<<endl;
	//  int i;
	//  cin>>i;      
    //}
    
    return coag;
  }
}

#define COEFFICIENT_REPARTITION_FILE_CLASS_PARTICLE_CXX
#endif
