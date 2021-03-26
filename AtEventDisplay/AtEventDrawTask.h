/**
 * @brief Event display task
 * @author JungWoo Lee (Korea Univ.)
 *         Adapted for AtTPCROOT by Yassid Ayyad (NSCL)
 */
#ifndef AtEVENTDRAWTASK_H
#define AtEVENTDRAWTASK_H

// FairRoot classes
#include "FairTask.h"
#include "FairLogger.h"

// ROOT classes
#include "TEvePointSet.h"
#include "TEveGeoShape.h"
#include "TEveBoxSet.h"
#include "TEveLine.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TPaletteAxis.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TH3.h"
#include "TGraph.h"
#include "TH2Poly.h"

#include "AtEventManager.h"
#include "AtRawEvent.h"
#include "AtEvent.h"
#include "AtProtoEvent.h"
#include "AtPatternEvent.h"
#include "AtTrackingEventAna.h"
#include "AtHoughSpaceLine.h"
#include "AtHoughSpaceCircle.h"
#include "AtHoughSpace.h"
#include "AtRansac.h"
#include "AtRansacMod.h"
#include "AtMlesacMod.h"
#include "AtLmedsMod.h"
#include "AtTrackFinderHC.h"
#include "AtHit.h"
#include "AtTpcMap.h"
#include "AtMap.h"
#include "AtTpcProtoMap.h"
#include "AtProtoQuadrant.h"
#include <fstream>

#ifndef __CINT__ // Boost
#include <boost/multi_array.hpp>
#endif //__CINT__

class AtEventDrawTask : public FairTask {
public:
   AtEventDrawTask();
   AtEventDrawTask(TString modes);

   virtual ~AtEventDrawTask();

   virtual InitStatus Init();
   virtual void Exec(Option_t *option);
   void Reset();

   // void Set2DPlotRange(Int_t uaIdx);
   void SetThreshold(Int_t val) { fThreshold = val; }
   void UnpackHoughSpace() { fUnpackHough = kTRUE; }
   void SetHitAttributes(Color_t, Size_t, Style_t);
   void Set3DHitStyleBar();
   void Set3DHitStyleBox();
   void SetSaveTextData();
   void SetLine(double t, std::vector<Double_t> p, double &x, double &y, double &z);
   void SetLine6(double t, std::vector<Double_t> p, double &x, double &y, double &z);
   // void SetHitClusterAttributes(Color_t, Size_t, Style_t);
   // void SetRiemannAttributes(Color_t, Size_t, Style_t);

   static void SelectPad(const char *rawevt);
   void DrawWave(Int_t PadNum);

   void SetGeoOption(Option_t *option) { fGeoOption = option; }

   void SetProtoMap(TString map) { fMap = map; }

   void SetMultiHit(Int_t hitMax);
   void SetAlgorithm(Int_t val) { fRANSACAlg = val; };

protected:
   virtual void DrawPadPlane();
   virtual void DrawPadWave();
   virtual void DrawPadAll();
   virtual void DrawQEvent();
   virtual void DrawRhoVariance();
   virtual void DrawHoughSpace();
   virtual void DrawHoughSpaceProto();
   virtual void DrawPhiReco();
   virtual void DrawMesh();
   virtual void Draw3DHist();
   virtual void DrawRad();
   virtual void DrawTheta();
   virtual void DrawThetaxPhi();
   virtual void DrawMC();
   virtual void DrawAux();

   AtMap *fAtMapPtr;
   void UpdateCvsPadPlane();
   void UpdateCvsPadWave();
   void UpdateCvsPadAll();
   void UpdateCvsQEvent();
   void UpdateCvsRhoVariance();
   void UpdateCvsHoughSpace();
   void UpdateCvsPhi();
   void UpdateCvsMesh();
   void UpdateCvs3DHist();
   void UpdateCvsRad();
   void UpdateCvsTheta();
   void UpdateCvsThetaxPhi();
   void UpdateCvsQuadrants();
   void UpdateCvsMC();
   void UpdateCvsAux();

   void ResetPadAll();
   void ResetPhiDistr();

   void DrawHitPoints();
   void DrawHSpace();
   void DrawProtoSpace();
   void DrawMeshSpace();
   // void DrawHitClusterPoints();
   // void DrawRiemannHits();

   EColor GetTrackColor(int i);

   Bool_t fIs2DPlotRange;
   Bool_t fUnpackHough;
   Bool_t fIsCircularHough;
   Bool_t fIsLinearHough;

   TClonesArray *fHitArray;
   TClonesArray *fRawEventArray;
   TClonesArray *fHoughSpaceArray;
   TClonesArray *fProtoEventArray;
   TClonesArray *fRansacArray;
   TClonesArray *fTrackFinderHCArray;
   TClonesArray *fTrackingEventAnaArray;
   TClonesArray *fPatternEventArray;

   AtHoughSpaceLine *fHoughSpaceLine_buff;
   AtHoughSpaceCircle *fHoughSpaceCircle_buff;
   AtRANSACN::AtRansac *fRansac;
   AtRansacMod *fRansacMod;
   AtMlesacMod *fMlesacMod;
   AtLmedsMod *fLmedsMod;
   AtTrackingEventAna *fTrackingEventAna;
   AtPATTERN::AtTrackFinderHC *fTrackFinderHC;

   AtEventManager *fEventManager;
   AtRawEvent *fRawevent;

   AtMap *fDetmap;

   Int_t fThreshold;
   Option_t *fGeoOption; // Chose Geometry of the detector: AtTPC (Default)-  Prototype
   TString fMap;

   TEvePointSet *fHitSet;
   TEvePointSet *fHitSetMin;

   TEvePointSet *fHitSetMC[5];    // For MC results
   TEvePointSet *fHitSetTFHC[10]; // for TrackFinderHC

   // TEveGeoShape* x;
   // std::vector<TEveGeoShape*> hitSphereArray;

   TEveBoxSet *fhitBoxSet;

   TPaletteAxis *fPadPlanePal;

   Color_t fHitColor;
   Size_t fHitSize;
   Style_t fHitStyle;

   /*TEvePointSet* fHitClusterSet;
   Color_t fHitClusterColor;
   Size_t  fHitClusterSize;
   Style_t fHitClusterStyle;*/

   /*vector<TEvePointSet*> fRiemannSetArray;
   Color_t fRiemannColor;
   Size_t  fRiemannSize;
   Style_t fRiemannStyle;*/

   TCanvas *fCvsPadPlane;
   TH2Poly *fPadPlane;
   TCanvas *fCvsPadWave;
   TH1I *fPadWave;
   TCanvas *fCvsPadAll;
   TH1I *fPadAll[300];
   TCanvas *fCvsQEvent;
   TH1D *fQEventHist;
   TH1D *fQEventHist_H;
   TCanvas *fCvsHoughSpace;
   TH2F *fHoughSpace;
   TCanvas *fCvsRhoVariance;
   TH1D *fRhoVariance;
   TCanvas *fCvsPhi;
   TH1D *fPhiDistr[5];
   TCanvas *fCvsMesh;
   TH1F *fMesh;
   TCanvas *fCvs3DHist;
   TH3F *f3DHist;
   TCanvas *fCvsRad;
   TH2F *fRadVSTb;
   TCanvas *fCvsTheta;
   TH2F *fTheta;
   TCanvas *fCvsThetaxPhi;
   TH2F *fThetaxPhi;
   TCanvas *fCvsQuadrant1;
   TH2F *fQuadrant1;
   TCanvas *fCvsQuadrant2;
   TH2F *fQuadrant2;
   TCanvas *fCvsQuadrant3;
   TH2F *fQuadrant3;
   TCanvas *fCvsQuadrant4;
   TH2F *fQuadrant4;
   TH1F *fAuxChannels[9];
   TCanvas *fCvsAux;

   TH2F *fThetaxPhi_Ini;
   TH2F *fThetaxPhi_Ini_RANSAC;

   TCanvas *fCvsMC_XY;
   TGraph *fMC_XY;
   TGraph *fMC_XY_exp;
   TGraph *fMC_XY_int;
   TGraph *fMC_XY_back;
   TCanvas *fCvsMC_Z;
   TGraph *fMC_ZX;
   TGraph *fMC_ZX_int;
   TGraph *fMC_ZX_back;
   TGraph *fMC_ZY;
   TGraph *fMC_ZY_int;
   TGraph *fMC_ZY_back;

   Int_t fNQuads;

   Int_t fMinZ;
   Int_t fMaxZ;
   Int_t fMinX;
   Int_t fMaxX;

   Int_t f3DHitStyle;
   Int_t fMultiHit;
   Bool_t fSaveTextData;
   Float_t f3DThreshold;
   Bool_t fIsRawData;
   Int_t fRANSACAlg;

   TF1 *fHoughLinearFit;
   TF1 *fRansacLinearFit;
   AtHit const *fIniHit;
   AtHit const *fIniHitRansac;

   // std::vector<TEveLine*> fLineArray;
   TEveLine *fLineArray[5];
   TEvePointSet *fVertex = nullptr;
   Int_t fLineNum;
   Int_t fTrackNum;
   // TEveLine* fLine;

   ClassDef(AtEventDrawTask, 1);
};

#endif
