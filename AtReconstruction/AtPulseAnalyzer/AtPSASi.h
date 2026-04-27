#ifndef AtPSASI_H
#define AtPSASI_H

#include "AtGenericTrace.h"
#include "AtPSA.h"

#include <Rtypes.h> // for Bool_t, THashConsistencyHolder, ClassDefOverride

#include <array>  // for array
#include <memory> // for make_unique, unique_ptr

class AtPad;
class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @brief Simple max or min finding PSA method for RCNP Si detectors.
 *
 *
 *
 */
class AtPSASi : public AtPSA {

private:
   Bool_t fIsTimeCorr{false};
   Bool_t fPositivePolarity{true};
   int low_bl_region[2];
   int hi_bl_region[2];
   int energy_integral[2];

public:
   AtPSASi();
   virtual HitVector AnalyzePad(AtPad *pad) override;
   std::unique_ptr<AtPSA> Clone() override { return std::make_unique<AtPSASi>(*this); }
   virtual HitVector AnalyzeGenTrace(AtGenericTrace *genTrace);

   void SetTimeCorrection(Bool_t value) { fIsTimeCorr = value; }
   void SetPositivePolarity(Bool_t value) { fPositivePolarity = value; }
   void SetLowBLRegion(int low, int hi);
   void SetHighBLRegion(int low, int hi);
   void SetEnergyIntegral(int low, int hi);

private:
   bool shouldSaveHit(double charge, double threshold, int tb);
   Double_t getTBCorr(std::array<Double_t, 512> &trace, int maxAdcIdx);

   ClassDefOverride(AtPSASi, 1)
};

#endif
