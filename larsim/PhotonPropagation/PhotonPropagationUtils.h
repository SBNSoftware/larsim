////////////////////////////////////////////////////////////////////////
// Description:
// Utility functions
////////////////////////////////////////////////////////////////////////
#ifndef PhotonPropagationUtils_H
#define PhotonPropagationUtils_H

#include <array>
#include <cmath>
#include <limits>
#include <vector>

/* #include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h" */

/* namespace fhicl { */
/*   class ParameterSet; */
/* } */

namespace phot
{
  /* class ScintTimeLAr : public ScintTime */
  /* { */
  /* public: */
  /*   explicit ScintTimeLAr(fhicl::ParameterSet const& pset); */
  /*   void GenScintTime(bool is_fast, CLHEP::HepRandomEngine& engine); */

  /* private: */
  /*   int           LogLevel; */

  /*   // parameters for the shape of argon scinitllation light time distribution */
  /*   double         SRTime;                        // PureLAr: rising time of slow LAr scinitllation; */
  /*   double         SDTime;                        // PureLAr: decay time of slow LAr scintillation; */
  /*   double         FRTime;                        // PureLAr: rising time of fast LAr scinitllation; */
  /*   double         FDTime;                        // PureLAr: decay time of fast LAr scintillation; */

  /*   // general functions */
  /*   double single_exp(double t, double tau2) const; */
  /*   double bi_exp(double t, double tau1, double tau2) const; */
  /* };
   */

  double fast_acos(double x);
  double interpolate(const std::vector<double>& xData,
                     const std::vector<double>& yData,
                     double x,
                     bool extrapolate,
                     size_t i = 0);
  double interpolate2(const std::vector<double>& xDistances,
                      const std::vector<double>& rDistances,
                      const std::vector<std::vector<std::vector<double>>>& parameters,
                      const double x,
                      const double r,
                      const size_t k);
  void interpolate3(std::array<double, 3>& inter,
                    const std::vector<double>& xData,
                    const std::vector<double>& yData1,
                    const std::vector<double>& yData2,
                    const std::vector<double>& yData3,
                    double x,
                    bool extrapolate);

    // implements relative method - do not use for comparing with zero
    // use this most of the time, tolerance needs to be meaningful in your context
    template <typename TReal>
    inline constexpr static bool
    isApproximatelyEqual(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
    {
        TReal diff = std::fabs(a - b);
        if (diff <= tolerance) return true;
        if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
        return false;
    }

    // supply tolerance that is meaningful in your context
    // for example, default tolerance may not work if you are comparing double with
    // float
    template <typename TReal>
    inline constexpr static bool
    isApproximatelyZero(TReal a, TReal tolerance = std::numeric_limits<TReal>::epsilon())
    {
        if (std::fabs(a) <= tolerance) return true;
        return false;
    }

    // use this when you want to be on safe side
    // for example, don't start rover unless signal is above 1
    template <typename TReal>
    inline constexpr static bool
    isDefinitelyLessThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
    {
        TReal diff = a - b;
        if (diff < tolerance) return true;
        if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
        return false;
    }

    template <typename TReal>
    inline constexpr static bool
    isDefinitelyGreaterThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
    {
        TReal diff = a - b;
        if (diff > tolerance) return true;
        if (diff > std::fmax(std::fabs(a), std::fabs(b)) * tolerance) return true;
        return false;
    }

}
#endif
