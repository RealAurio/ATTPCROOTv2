#include "AtTabMAGNEX.h"

#include "AtTabMain.h"
#include "AtHitClusterEvent.h"
#include "AtViewerManager.h"

#include <FairLogger.h>

#include <TEveManager.h>

namespace DataHandling {
class AtSubject;
}

ClassImp(AtTabMAGNEX);

AtTabMAGNEX::AtTabMAGNEX() : AtTabMain()
{
   fHitClusterEventBranch = &AtViewerManager::Instance()->GetHitClusterEventBranch();
   fHitClusterEventBranch->Attach(this);
}


AtTabMAGNEX::~AtTabMAGNEX()
{
   fHitClusterEventBranch->Detach(this);
   AtTabMain::~AtTabMain();
}

void AtTabMAGNEX::InitTab()
{
   AtTabMain::InitTab();

   gEve->AddEvent(fEveHitClusterEvent.get());
   fTabInfo->AddAugment(std::make_unique<AtTabInfoFairRoot<AtHitClusterEvent>>(*fHitClusterEventBranch));
}

void AtTabMAGNEX::Update(DataHandling::AtSubject *sub)
{
   if (sub == fHitClusterEventBranch || sub == fEntry)
      UpdateHitClusterEventElements();

   AtTabMain::Update(sub);
}

void AtTabMAGNEX::UpdateRenderState()
{
   fEveHitClusterEvent->SetRnrState(true);
   AtTabMain::UpdateRenderState();
}

void AtTabMAGNEX::UpdateHitClusterEventElements()
{
   if (fEveHitClusterEvent == nullptr)
      return;

   auto fHitClusterEvent = GetFairRootInfo<AtHitClusterEvent>();
   if (fHitClusterEvent == nullptr) {
      LOG(debug) << "Cannot update AtHitClusterEvent elements: no event available";
      return;
   }
   LOG(info) << "Hello there.";

}
