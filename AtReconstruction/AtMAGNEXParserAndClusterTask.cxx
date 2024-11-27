#include "AtMAGNEXParserAndClusterTask.h"

#include <FairTask.h>
#include <FairLogger.h>
#include <FairRootManager.h>

#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <Math/Point2D.h>
#include <Math/Point3D.h>

#include <memory>
#include <map>
#include <fstream>
#include <utility>
#include <string>
#include <sstream>
#include <vector>

#include "AtHitClusterEvent.h"
#include "AtHitClusterFull.h"

using XYPoint = ROOT::Math::XYPoint;
using XYZPoint = ROOT::Math::XYZPoint;

AtMAGNEXParserAndClusterTask::AtMAGNEXParserAndClusterTask(TString inputFileName, TString channelToColMapFileName, TString detParFileName, TString outputBranchName)
   : fInputFileName(std::move(inputFileName)), fOutputBranchName(std::move(outputBranchName)),
     fChannelToColMapFileName(std::move(channelToColMapFileName)), fDetectorParametersFileName(std::move(detParFileName)),
     fEventArray("AtHitClusterEvent", 1)
{
}

InitStatus AtMAGNEXParserAndClusterTask::Init()
{
   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not find RootManager!";
      return kERROR;
   }
   ioMan->Register(fOutputBranchName, "MAGNEX", &fEventArray, fIsPersistence);

   fInputFile = std::make_unique<TFile>(fInputFileName, "READ");
   if (fInputFile->IsZombie()) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not open ROOT file " << fInputFileName << "!";
      return kERROR;
   } else {
      LOG(info) << "ROOT file " << fInputFileName << " opened successfully.";
   }

   fInputTree = (TTree *)fInputFile->Get("Data_R");
   LOG(info) << "Total number of entries in the tree: " << fInputTree->GetEntries() << ".";
   fInputTree->SetBranchAddress("Board", &Board);
   fInputTree->SetBranchAddress("Channel", &Channel);
   fInputTree->SetBranchAddress("FineTSInt", &FTS);
   fInputTree->SetBranchAddress("CoarseTSInt", &CTS);
   fInputTree->SetBranchAddress("Timestamp", &Timestamp);
   fInputTree->SetBranchAddress("Charge", &Charge);
   fInputTree->SetBranchAddress("Flags", &Flags);
   fInputTree->SetBranchAddress("Pads", &Pad);
   fInputTree->SetBranchAddress("Charge_cal", &Charge_cal);
   fInputTree->SetBranchAddress("Row", &Row);
   fInputTree->SetBranchAddress("Section", &Section);

   fChannelToColMapFile = std::make_unique<std::ifstream>(fChannelToColMapFileName);
   fChannelToColMap = std::make_unique<std::map<Int_t, Int_t>>();
   if (!fChannelToColMapFile->is_open()) {
      LOG(error) << "AtMAGNEXParserAndClusterTask::Init Error : Could not open the channel to column (strip) file " << fChannelToColMapFileName << "!";
      return kERROR;
   } else {
      LOG(info) << "Channel to column (strip) file " << fChannelToColMapFileName << " opened successfully.";

      std::string line;
      for (Int_t i = 0; i < 8; ++i) {
         std::getline(*fChannelToColMapFile, line);
         LOG(debug) << line;
      }

      while (!fChannelToColMapFile->eof()) {
         Int_t channel, column;
         std::getline(*fChannelToColMapFile, line);
         std::istringstream buffer(line);
         buffer >> channel >> column;
         LOG(debug) << "Channel: " << channel << " Column: " << column << ".";
         fChannelToColMap->insert(std::make_pair(channel, column));
      }
   }

   fMap = new AtMAGNEXMap(fDetectorParametersFileName);
   fMap->GeneratePadPlane();

   return kSUCCESS;
}

void AtMAGNEXParserAndClusterTask::Exec(Option_t *opt)
{
   fInputTree->GetEntry(fEntryNum);
   ULong64_t timeinit = Timestamp;

   auto *event = dynamic_cast<AtHitClusterEvent *>(fEventArray.ConstructedAt(0, "C"));
   event->SetEventID(fEventNum++);
   event->SetTimestamp(timeinit);

   std::vector<std::vector<ULong64_t>> entriesByRow;
   for (Int_t i = 0; i < fMap->GetRowNum(); ++i) {
      std::vector<ULong64_t> entriesVector{};
      entriesByRow.push_back(entriesVector);
   }

   while (fEntryNum < fInputTree->GetEntries()) {
      fInputTree->GetEntry(fEntryNum);
      if (Timestamp - timeinit > fWindowSize)
         break;
      Int_t Col = fChannelToColMap->find(Channel)->second;
      if (Col < fMap->GetColNum())
         entriesByRow[Row].push_back(fEntryNum);
      ++fEntryNum;
   }

   Int_t hitNum{0}, clusterNum{0};
   for (Int_t i = 0; i < fMap->GetRowNum(); ++i) {
      if (entriesByRow[i].size() == 0)
         continue;

      AtHitClusterFull *hitCluster = new AtHitClusterFull();
      hitCluster->SetClusterID(clusterNum++);

      for (ULong64_t entry : entriesByRow[i]) {
         fInputTree->GetEntry(entry);

         Int_t Col = fChannelToColMap->find(Channel)->second;

         XYPoint point = fMap->CalcPadCenter(fMap->PadID(Col, Row));
         Double_t x = point.X();
         Double_t y = FTS;               // Still needs calibration.
         Double_t z = point.Y();

         AtHit *hit = new AtHit(hitNum++, fMap->PadID(Col, Row), XYZPoint(x, 0, z), Charge);
         hit->SetTimeStamp(FTS);

         hitCluster->AddHit(*hit);
      }
      event->AddHitCluster(hitCluster);
   }

   LOG(info) << "Event " << fEventNum << " with timestamp " << timeinit << " has " << hitNum << " hits in " << clusterNum << " clusters.";

}

ClassImp(AtMAGNEXParserAndClusterTask)
