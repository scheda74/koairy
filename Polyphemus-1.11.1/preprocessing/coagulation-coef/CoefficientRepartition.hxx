
#ifndef COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_HXX

#include "CoefficientRepartitionHeader.hxx"

#include "talos/Talos.hxx"

#include "seldon/Seldon.hxx"

#include "Ops.hxx"

#include "Error.cxx"
#include "ClassCoefficientRepartitionBase.cxx"
#include "ClassParticle.cxx"
#include "ClassGeneralSection.cxx"
#include "ClassCoefficientRepartition.cxx"

namespace CoefficientRepartition
{
  // Number of species.
  int ClassCoefficientRepartitionBase::Ns_ = 1;

  int ClassCoefficientRepartitionBase::rank_ = 0;
  int ClassCoefficientRepartitionBase::Nrank_ = 1;

  // Density in µg/µm3
  double ClassCoefficientRepartitionBase::density_ = double(1e-6);

  // Configuration file.
  string ClassCoefficientRepartitionBase::configuration_file_ = "";

  // Number of size sections.
  int ClassGeneralSection::Nb_ = 0;

  // Number of composition sections.
  int ClassGeneralSection::Nc_ = 1;

  // Bound section diameters.
  Vector2T ClassGeneralSection::diameter_;

  // Bound section masses.
  Vector2T ClassGeneralSection::mass_;

  // Vectors of composition index and fractions.
  Vector3I ClassGeneralSection::fraction_index_;
  Vector3T ClassGeneralSection::fraction_;
}

#define COEFFICIENT_REPARTITION_FILE_COEFFICIENT_REPARTITION_HXX
#endif
