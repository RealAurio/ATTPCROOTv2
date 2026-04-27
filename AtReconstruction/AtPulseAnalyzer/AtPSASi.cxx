#include "AtPSASi.h"

#include "AtHit.h"
#include "AtPad.h"

#include <FairLogger.h>

#include <Math/Point3D.h>    // for PositionVector3D
#include <Math/Point3Dfwd.h> // for XYZPoint

#include <algorithm>
#include <array>    // for array
#include <iterator> // for distance
#include <memory>   // for unique_ptr, make_unique
#include <numeric>
#include <utility> // for pair

// #ifdef _OPENMP
// #include <omp.h>
// #endif
using XYZPoint = ROOT::Math::XYZPoint;

AtPSASi::AtPSASi() : AtPSA() { 
  low_bl_region[0] = -50;
  low_bl_region[1] = -30;
  hi_bl_region[0] = 80;
  hi_bl_region[1] = 100;
  energy_integral[0] = -10;
  energy_integral[1] = 15;
}

void AtPSASi::SetLowBLRegion(int low, int hi) {
  low_bl_region[0] = low;
  low_bl_region[1] = hi;
}

void AtPSASi::SetHighBLRegion(int low, int hi) {
  hi_bl_region[0] = low;
  hi_bl_region[1] = hi;
}

void AtPSASi::SetEnergyIntegral(int low, int hi) {
  energy_integral[0] = low;
  energy_integral[1] = hi;
}

AtPSASi::HitVector AtPSASi::AnalyzePad(AtPad *pad)
{
   XYZPoint pos(0, 0, 0);

   std::array<Double_t, 512> floatADC = pad->GetADC();

   // Get baseline value.
   double baseline{};
   Double_t *maxAdcIt{nullptr};
   Double_t charge{-999};
   Int_t maxAdcIdx;
   if (fPositivePolarity) {
      maxAdcIt = std::max_element(floatADC.begin() + 20, floatADC.end() - 12);
      maxAdcIdx = std::distance(floatADC.begin(), maxAdcIt);
   } else {
      maxAdcIt = std::min_element(floatADC.begin() + 20, floatADC.end() - 12);
      maxAdcIdx = std::distance(floatADC.begin(), maxAdcIt);
      (*maxAdcIt) = -(*maxAdcIt);
   }

   float count = 0;
   for (int i = std::max(20, maxAdcIdx + low_bl_region[0]); i < std::min(maxAdcIdx + low_bl_region[1], 500);
        i++) { // more baseline added
      baseline += floatADC[i];
      count += 1;
   }
   for (int i = std::max(20, maxAdcIdx + hi_bl_region[0]); i < std::min(maxAdcIdx + hi_bl_region[1], 500);
        i++) { // more baseline added
      baseline += floatADC[i];
      count += 1;
   }
   baseline /= count;
   LOG(debug) << "Baseline calculated: " << baseline;
   charge = *maxAdcIt - baseline;

   if (!shouldSaveHit(charge, fThreshold, maxAdcIdx))
      return {};

   LOG(debug) << "========== float ADC array max value ========== " << *maxAdcIt;
   LOG(debug) << "========== float ADC array max index: " << maxAdcIdx << " ==========";

   // Calculation of the mean value of the peak time by interpolating the pulse
   Double_t timemax = 0.5 * (floatADC[maxAdcIdx - 1] - floatADC[maxAdcIdx + 1]) /
                      (floatADC[maxAdcIdx - 1] + floatADC[maxAdcIdx + 1] - 2 * floatADC[maxAdcIdx]);
   Double_t TBCorr = getTBCorr(floatADC, maxAdcIdx);
   Double_t QHitTot = std::abs(std::accumulate(floatADC.begin() + maxAdcIdx + energy_integral[0],
                                               floatADC.begin() + maxAdcIdx + energy_integral[1], 0) /
                                  (energy_integral[1] - energy_integral[0]) -
                               baseline);

   auto hit = std::make_unique<AtHit>(0, pos, charge);

   hit->SetTimeStamp(maxAdcIdx);
   hit->SetTimeStampCorr(TBCorr);
   hit->SetTimeStampCorrInter(timemax);
   hit->SetTraceIntegral(QHitTot);

   HitVector ret;
   ret.push_back(std::move(hit));
   return ret;
}
bool AtPSASi::shouldSaveHit(double charge, double threshold, int tb)
{
   bool ret = true;
   LOG(debug) << charge << "   " << threshold << "  " << tb;
   if (threshold > 0 && charge < threshold) {
      ret = false;
      LOG(debug) << "Invalid threshold with charge: " << charge << " and threshold: " << threshold;
   }
   if ((tb < 20 || tb > 500)) {
      ret = false;
      LOG(debug) << "Peak is outside valid time window (20,500) TBs.";
   }

   return ret;
}

Double_t AtPSASi::getTBCorr(AtPad::trace &adc, int maxAdcIdx)
{
   if (maxAdcIdx < 11)
      return 0;

   Double_t qTot = 0;
   Double_t tbAvg = 0;
   for (Int_t i = 0; i < 11; i++)
      if (adc[maxAdcIdx - i + 10] > 0 && adc[maxAdcIdx - i + 10] < 4000) {
         auto tb = maxAdcIdx - i + 5;
         qTot += adc[tb];
         tbAvg += adc[tb] / qTot * (tb - tbAvg);
      }
   return tbAvg;
}

AtPSASi::HitVector AtPSASi::AnalyzeGenTrace(AtGenericTrace *genTrace)
{
   XYZPoint pos(0, 0, 0);

   std::array<Double_t, 256> floatADC;
   std::vector<Double_t> floatADCVector = genTrace->GetADC();
   // FOR GAGG
   //for (int i = 0; i < 256; i++)
   //   floatADC[i] = floatADCVector[i];
   // std::cout << "GAGG ADC entries: " << floatADCVector.size() <<  std::endl;
   if (floatADCVector.size() >= 256) {
      for (int i = 0; i < 256; i++) {
         floatADC[i] = floatADCVector[i];
      }
   } else {
      LOG(error) << "There are not 256 ADC values in the GAGG trace. Skipping!";
      return {};
   }

   Double_t *maxAdcIt{nullptr};
   Double_t charge{-999};
   Double_t baseline{0};
   if (fPositivePolarity) {
      maxAdcIt = std::max_element(floatADC.begin() + 20, floatADC.end() - 12);
      charge = *maxAdcIt;
   } else {
      maxAdcIt = std::min_element(floatADC.begin() + 20, floatADC.end() - 12);
      charge = -(*maxAdcIt);
   }
   Int_t maxAdcIdx = std::distance(
      floatADC.begin(),
      maxAdcIt); // To-Do: might need to utilize similar trace intergal calculation as used in Si. Unsure yet.

   float count = 0;
   for (int i = std::max(20, maxAdcIdx + low_bl_region[0]); i < std::min(maxAdcIdx + low_bl_region[1], (int)(floatADCVector.size()-12));
        i++) { // more baseline added
      baseline += floatADC[i];
      count += 1;
   }
   for (int i = std::max(20, maxAdcIdx + hi_bl_region[0]); i < std::min(maxAdcIdx + hi_bl_region[1], (int)(floatADCVector.size()-12));
        i++) { // more baseline added
      baseline += floatADC[i];
      count += 1;
   }
   baseline /= count;
   LOG(debug) << "Baseline calculated: " << baseline;
   charge = *maxAdcIt - baseline;

   //std::cout << " Max ADC = " << *maxAdcIt << std::endl;

   //std::cout << " Diff ADC = " << *maxAdcIt - baseline << std::endl;

   //if (abs(*maxAdcIt - baseline) > 100) {
   //   for (int i=0; i<256; ++i) {
   //     std::cout << i << "  " << floatADC[i] << std::endl;
   //   }
   //}

   if (!shouldSaveHit(charge, fThreshold, maxAdcIdx)) {
      LOG(debug) << "GAGG trace did not pass threshold.";
      return {};
   }

   // Calculation of the mean value of the peak time by interpolating the pulse
   Double_t timemax = 0.5 * (floatADC[maxAdcIdx - 1] - floatADC[maxAdcIdx + 1]) /
                      (floatADC[maxAdcIdx - 1] + floatADC[maxAdcIdx + 1] - 2 * floatADC[maxAdcIdx]);
   Double_t QHitTot = std::abs(std::accumulate(floatADC.begin() + maxAdcIdx + energy_integral[0],
                                               floatADC.begin() + maxAdcIdx + energy_integral[1], 0) /
                                  (energy_integral[1] - energy_integral[0]) -
                               baseline);

   auto hit = std::make_unique<AtHit>(0, pos, charge);

   hit->SetTimeStamp(maxAdcIdx);
   hit->SetTimeStampCorrInter(timemax);
   hit->SetTraceIntegral(QHitTot);

   HitVector ret;
   ret.push_back(std::move(hit));
   return ret;
}

ClassImp(AtPSASi);
