#include "AtEvent.h"
#include "AtHit.h"
#include "AtMCPoint.h"
#include "AtPad.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRun.h>
#include <FairRunAna.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMaterialInterface.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSpectrum.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreePlayer.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TrackPoint.h>

#include "AbsFitterInfo.h"
#include "AbsKalmanFitter.h"
#include "ConstField.h"
#include "DAF.h"
#include "EventDisplay.h"
#include "Exception.h"
#include "FieldManager.h"
#include "FitStatus.h"
#include "KalmanFitStatus.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitterRefTrack.h"
#include "MaterialEffects.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementFactory.h"
#include "MeasurementOnPlane.h"
#include "MeasurementProducer.h"

#include <ios>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <vector>

// ROOT
#include <TApplication.h>
#include <TArrayD.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRotation.h>
#include <TVectorD.h>

#include "Fit/Fitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/EulerAngles.h"
#include "Math/GenVector/Quaternion.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationX.h"
#include "Math/GenVector/RotationY.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/Minimizer.h"

struct trackSegment {
   Double_t eLoss;
   TVector3 iniPos;
   TVector3 deltaMom;
   TVector3 deltaPos;
   Double_t theta;
   Double_t phi;
   UInt_t id;

   friend std::ostream &operator<<(std::ostream &os, const trackSegment &ts);
};

std::ostream &operator<<(std::ostream &os, const trackSegment &ts)
{
   os << "\n";
   os << " Track segment : " << ts.id << " - Momentum:  " << ts.deltaMom.X() << " - " << ts.deltaMom.Y() << " - "
      << ts.deltaMom.Z() << " - Energy Loss : " << ts.eLoss << "\n";
   os << " =============   - Position :  " << ts.iniPos.X() << " - " << ts.iniPos.Y() << " - " << ts.iniPos.Z()
      << " . Mag Dir : " << ts.iniPos.Mag() << "\n";
   os << " =============   - Position direction :  " << ts.deltaPos.X() << " - " << ts.deltaPos.Y() << " - "
      << ts.deltaPos.Z() << " . Mag Dir : " << ts.deltaPos.Mag() << "\n";
   os << " =============   - Theta    :  " << ts.theta * TMath::RadToDeg() << " - Phi : " << ts.phi * TMath::RadToDeg()
      << "\n";
   return os;
}

struct firstOrbit {
   Double_t POCA;
   Double_t Z;
   Double_t phi;
   Double_t length;
   Double_t eLoss;
};

std::tuple<Double_t, Double_t>
GetMomFromBrho(Double_t A, Double_t Z, Double_t brho); ///< Returns momentum (in GeV) from Brho assuming M (amu) and Z;
double
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);
Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
Double_t GetNPeaksHRS(std::vector<Int_t> *timeMax, std::vector<Float_t> *adcMax, double *adc_test);
double GetMaximum(double *adc);
void ClusterizeSmooth3D(AtTrack &track, Float_t distance, Float_t radius);
void Clusterize3D(AtTrack &track, Float_t distance, Float_t radius);
void Clusterize(AtTrack &track);
firstOrbit GetFirstOrbit(genfit::Track *track, genfit::AbsTrackRep *rep, TVector3 vertex);
void ConstructTrack(const genfit::StateOnPlane *prevState, const genfit::StateOnPlane *state,
                    const genfit::AbsTrackRep *rep, std::vector<TVector3> &track, std::vector<trackSegment> &segments);
