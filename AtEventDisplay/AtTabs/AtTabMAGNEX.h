#ifndef ATTABMAGNEX_H
#define ATTABMAGNEX_H

#include "AtTabMain.h"
#include "AtViewerManagerSubject.h"
#include "AtHitClusterFull.h"

#include <Rtypes.h>

class TBuffer;
class TClass;
class TMemberInspector;
namespace DataHandling {
class AtSubject;
}

class AtTabMAGNEX : public AtTabMain {
protected:
   TEveEventManagerPtr fEveHitClusterEvent{std::make_unique<TEveEventManager>("AtHitClusterEvent")};
   std::vector<TEvePointSetPtr> fHitClusterSets;

   DataHandling::AtBranch *fHitClusterEventBranch;

public:
   AtTabMAGNEX();
   ~AtTabMAGNEX();

   void InitTab() override;
   void Update(DataHandling::AtSubject *sub) override;

protected:
   void UpdateRenderState() override;

   void SetPointsFromHitCluster(TEvePointSet &hitSet, const AtHitClusterFull &hitCluster);

private:
   void UpdateHitClusterEventElements();

   void ExpandNumHitClusters(Int_t num);

   ClassDefOverride(AtTabMAGNEX, 1);
};

#endif
