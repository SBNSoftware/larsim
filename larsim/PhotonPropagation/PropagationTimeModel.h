#ifndef PROPAGATIONTIMEMODEL_H
#define PROPAGATIONTIMEMODEL_H

// PropagationTimeModel
//  - parameterized fast optical simulation of photon propagation times
//  - contains functions to calculate the propagation times of direct and
//  reflected photons incident each photo-detector, along with the necessary
//  ultility functions
//  - full description of model: Eur. Phys. J. C 81, 349 (2021)

// Nov 2021 by P. Green

// LArSoft libraries
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larsim/IonizationScintillation/ISTPC.h"

// fhicl
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"

namespace CLHEP {
class HepRandomEngine;
}

// ROOT headers
#include "TF1.h"

#include <array>
#include <vector>

namespace phot {

class PropagationTimeModel {

public:
  // constructor
  PropagationTimeModel(fhicl::ParameterSet VUVTimingParams,
                       fhicl::ParameterSet VISTimingParams,
                       CLHEP::HepRandomEngine &scintTimeEngine,
                       bool doReflectedLight = false,
                       bool GeoPropTimeOnly = false);

  // propagation time
  void propagationTime(std::vector<double> &arrival_time_dist,
                       geo::Point_t const &x0, const size_t OpChannel,
                       bool Reflected = false);

private:
  void generateVUVParams(fhicl::ParameterSet &VUVTimingParams,
                         CLHEP::HepRandomEngine &scintTimeEngine);

  // direct / VUV light
  void getVUVTimes(std::vector<double> &arrivalTimes,
                   const double distance_in_cm, const size_t angle_bin);

  void getVUVTimesGeo(std::vector<double> &arrivalTimes,
                      const double distance_in_cm) const;

  // reflected / visible light
  void getVISTimes(std::vector<double> &arrivalTimes,
                   const geo::Point_t &scintPoint, const geo::Point_t &opDetPoint);

  // utility functions
  static double finter_d(const double *x, const double *par);

  static double model_close(const double *x, const double *par);

  static double model_far(const double *x, const double *par);

  // configuration
  const bool fGeoPropTimeOnly;

  // random numbers
  CLHEP::RandFlat fUniformGen;

  // geometry properties
  double fplane_depth;
  geo::Point_t fcathode_centre;

  // photodetector geometry properties
  std::vector<geo::Point_t> fOpDetCenter;
  std::vector<int> fOpDetOrientation;

  // For VUV propagation time parametrization
  double fstep_size, fvuv_vgroup_mean, fvuv_vgroup_max, fmin_d,
      finflexion_point_distance, fangle_bin_timing_vuv;
  // vector containing generated VUV timing parameterisations
  std::vector<std::vector<CLHEP::RandGeneral>> fVUVTimingGen;
  // vector containing min and max range VUV timing parameterisations are
  // sampled to
  std::vector<std::vector<double>> fVUV_max;
  std::vector<std::vector<double>> fVUV_min;

  // For VIS propagation time parameterisation
  double fvis_vmean, fangle_bin_timing_vis;
  std::vector<double> fdistances_refl;
  std::vector<double> fradial_distances_refl;
  std::vector<std::vector<std::vector<double>>> fcut_off_pars;
  std::vector<std::vector<std::vector<double>>> ftau_pars;
};

} // namespace phot

#endif
