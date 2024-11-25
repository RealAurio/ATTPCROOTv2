#ifndef ATMAGNEXPARSERANDCLUSTERTASK_H
#define ATMAGNEXPARSERANDCLUSTERTASK_H

#include <FairTask.h>

#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <memory>
#include <map>
#include <fstream>

#include "AtMAGNEXMap.h"

class TBuffer;
class TClass;
class TMemberInspector;

class AtMAGNEXParserAndClusterTask : public FairTask {
protected:
   Bool_t fIsPersistence{false};
   TString fInputFileName;
   TString fOutputBranchName;

   TString fChannelToColMapFileName;
   std::unique_ptr<std::ifstream> fChannelToColMapFile;
   std::unique_ptr<std::map<Int_t, Int_t>> fChannelToColMap;

   TString fDetectorParametersFileName;
   AtMAGNEXMap *fMap;

   std::unique_ptr<TFile> fInputFile{nullptr};
   TTree *fInputTree{nullptr};
   TClonesArray fEventArray;

   Int_t fEventNum{0};
   ULong64_t fEntryNum{0};
   ULong64_t fWindowSize{2000000};

   UShort_t Board;
   UShort_t Channel;
   UShort_t FTS;
   ULong64_t CTS;
   ULong64_t Timestamp;
   UShort_t Charge;
   UInt_t Flags;
   UShort_t Pad;
   Double_t Charge_cal;
   UShort_t Row;
   UShort_t Section;

public:
   AtMAGNEXParserAndClusterTask(TString inputFileName, TString channelToColMapFileName, TString detParFileName, TString outputBranchName = "AtHitClusterEventH");

   void SetPersistence(bool value) { fIsPersistence = value; }
   void SetWindowSize(ULong64_t value) { fWindowSize = value; }

   virtual InitStatus Init() override;
   virtual void Exec(Option_t *opt) override;

   ClassDefOverride(AtMAGNEXParserAndClusterTask, 1);
};

#endif
