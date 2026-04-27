/*********************************************************************
 *   ATSI Mapping Class	AtSiMap.h			             *
 *   Author: Y. Ayyad            				     *
 *   Log: 13-02-2015 17:16 JST					     *
 *								     *
 *********************************************************************/

#ifndef ATGAGGMAP_H
#define ATGAGGMAP_H
#include "AtMap.h"

#include <Math/Point2Dfwd.h>
#include <Rtypes.h>
#include <TXMLDocument.h>
#include <TXMLNode.h>

class TBuffer;
class TClass;
class TMemberInspector;

class AtGAGGMap : public AtMap {

public:
   AtGAGGMap();
   ~AtGAGGMap();

   virtual void Dump() override;
   virtual void GeneratePadPlane() override;
   virtual ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;
   virtual Int_t BinToPad(Int_t binval) override { return binval - 1; };
   virtual void ParseAtTPCMap(TXMLNode *node) override;

   Int_t InhibitStrips(TString stripsFilePath);

   ClassDefOverride(AtGAGGMap, 1);

protected:
   Int_t fill_coord(int pindex, float padxoff, float padyoff, float triside, float fort);
};

#endif
