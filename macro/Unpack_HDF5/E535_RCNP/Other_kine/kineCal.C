TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);

void kineCal()
{
   FairRunAna *run = new FairRunAna(); // Forcing a dummy run
   //Bypassing cutfiles edit It just creates a giant cut that allows everything through
   auto makeOpenCut = []() {
      auto *c = new TCutG();
      c->SetVarX("x"); 
      c->SetVarY("y");
      c->SetPoint(0, -1e9, -1e9);
      c->SetPoint(1,  1e9, -1e9);
      c->SetPoint(2,  1e9,  1e9);
      c->SetPoint(3, -1e9,  1e9);
      c->SetPoint(4, -1e9, -1e9);
      return c;
   };
   //   TString outfname="./canvas_kine.root";
   TString outfname="Cal_hists_veto_v2.root"; //HERE AMBER change name to match files completed
   TFile *outfile=new TFile(outfname,"recreate");

   // AtMap to check if a hit belong to a big pad or small pad.
   TString scriptfile = "rcnp_map_size.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   AtTpcMap *map = new AtTpcMap();
   map->ParseXMLMap(mapDir.Data());
   map->GeneratePadPlane();

   // Punch through filter.
   double punchThroughThreshold = 20;
   AtTools::AtPunchThroughChecker punchThroughChecker = AtTools::AtPunchThroughChecker();
   punchThroughChecker.SetDistanceThreshold(punchThroughThreshold);

   // ELoss model for kinetic energy estimations.
   double density = 1.3378e-3; // 470Torr
   std::vector<std::tuple<int, int, int>> materialComponents;
   materialComponents.push_back(std::make_tuple(12, 6, 3));
   materialComponents.push_back(std::make_tuple(2, 1, 8));

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_p = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_p");
   eLossModelC3D8_p->SetMaterial(materialComponents);
   eLossModelC3D8_p->SetProjectile(1, 1, 1.007825031898);
   eLossModelC3D8_p->SetPDGCode("1000010010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_d = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_d");
   eLossModelC3D8_d->SetMaterial(materialComponents);
   eLossModelC3D8_d->SetProjectile(2, 1, 2.0135532);
   eLossModelC3D8_d->SetPDGCode("1000020010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_3He = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_3He");
   eLossModelC3D8_3He->SetMaterial(materialComponents);
   eLossModelC3D8_3He->SetProjectile(3, 2, 3.01602932197);
   eLossModelC3D8_3He->SetPDGCode("1000030020");

   //Bypass cut files edit
   TCutG *cutPIDproton  = makeOpenCut();
   TCutG *cutPIDdeuteron = makeOpenCut();
   TCutG *cutPID3He     = makeOpenCut();

   TCutG *cutSiB   = makeOpenCut();
   TCutG *cutSiC   = makeOpenCut();
   TCutG *cutSi13B = makeOpenCut();
   TCutG *cutSiBe  = makeOpenCut();

   TCutG *cutGaggB  = makeOpenCut();
   TCutG *cutGaggBe = makeOpenCut();
   TCutG *cutGaggLi = makeOpenCut();
   TCutG *cutGaggHe = makeOpenCut();

   
   // Kinematic curve
   TGraph* kinecurve_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");
   kinecurve_3HeGS->SetLineWidth(2);
   kinecurve_3HeGS->SetLineColor(kRed);
   TGraph* ang_lab_cm_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");// 5 deg pitch in theta_cm
   ang_lab_cm_3HeGS->SetMarkerStyle(8);
   ang_lab_cm_3HeGS->SetMarkerSize(1);
   TLegend *legend = new TLegend(0.6,0.2,0.85,0.5);
   legend->AddEntry(kinecurve_3HeGS,"3He G.S.","l");
   legend->AddEntry(ang_lab_cm_3HeGS,"5 deg pitch in #theta_{cm}","p");
   legend->SetFillColor(0);

   outfile->cd();
   
   // Histogram definitions.
   TH2F *histRangeVThetaLAB = new TH2F("histRangeVThetaLAB", "histRangeVThetaLAB", 180, 0, 180, 1030, 0, 1030);
   TH2F *histTrackThetaLABRange3He = new TH2F("histTrackThetaLABRange3He", "histTrackThetaLABRange3He", 1030, 0, 1030, 180, 0, 180);
   TH2F *histTrackThetaLABRange3HeSiGaggBe = new TH2F("histTrackThetaLABRange3HeSiGaggBe", "histTrackThetaLABRange3HeSiGaggBe", 506, 0, 1030, 90, 0, 180);
   TH2F *histEstimatedKinEVThetaLABTotal = new TH2F("histEstimatedKinEVThetaLABTotal", "histEstimatedKinEVThetaLABTotal", 360, 0, 180, 500, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3He = new TH2F("histEstimatedKinEVThetaLAB3He", "histEstimatedKinEVThetaLAB3He", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3HeSiBe = new TH2F("histEstimatedKinEVThetaLAB3HeSiBe", "histEstimatedKinEVThetaLAB3HeSiBe", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBe = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBe", "histEstimatedKinEVThetaLAB3HeSiGaggBe", 360, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBetrack1 = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBetrack1", "histEstimatedKinEVThetaLAB3HeSiGaggBetrack1", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBetrack2 = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBetrack2", "histEstimatedKinEVThetaLAB3HeSiGaggBetrack2", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLABkinett = new TH2F("histEstimatedKinEVThetaLABkinett", "histEstimatedKinEVThetaLABkinett", 180, 0, 180, 300, 0, 30);
   TH2F *histThetaLABThetaLAB = new TH2F("histThetaLABThetaLAB", "histThetaLABThetaLAB", 360, 0, 180, 360, 0, 180);
   TH2F *histThetaLABThetaLABKine = new TH2F("histThetaLABThetaLABKine", "histThetaLABThetaLABKine", 500, 0, 50, 150, 0, 15);
   TH2F *histThetaLABThetaLAB3HeSiGaggBe = new TH2F("histThetaLABThetaLAB3HeSiGaggBe", "histThetaLABThetaLAB3HeSiGaggBe", 360, 0, 180, 360, 0, 180);
   TH2F *histEstimatedKinEVThetaLAB2H  = new TH2F("histEstimatedKinEVThetaLAB2H", "histEstimatedKinEVThetaLAB2H", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB1H  = new TH2F("histEstimatedKinEVThetaLAB1H", "histEstimatedKinEVThetaLAB1H", 180, 0, 180, 250, 0, 20);      
   TH2F *histEstimatedKinEVThetaLABSiB  = new TH2F("histEstimatedKinEVThetaLABSiB", "histEstimatedKinEVThetaLABSiB", 180, 0, 180, 250, 0, 20);      
   TH2F *histEstimatedKinEVThetaLABSiC  = new TH2F("histEstimatedKinEVThetaLABSiC", "histEstimatedKinEVThetaLABSiC", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSi13B  = new TH2F("histEstimatedKinEVThetaLABSi13B", "histEstimatedKinEVThetaLABSi13B", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSiBe  = new TH2F("histEstimatedKinEVThetaLABSiBe", "histEstimatedKinEVThetaLABSiBe", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSiGaggBe  = new TH2F("histEstimatedKinEVThetaLABSiGaggBe", "histEstimatedKinEVThetaLABSiGaggBe", 180, 0, 180, 250, 0, 20);
   TH2F *histdEdxVTotalRange = new TH2F("histdEdxVTotalRange", "histdEdxVTotalRange", 500, 0, 1030, 1600, 0, 4000);
   TH2F *histESmallVTotalRange = new TH2F("histESmallVTotalRange", "histESmallVTotalRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *histEBigVBigRange = new TH2F("histEBigVBigRange", "histEBigVBigRange", 500, 0, 1030, 1600, 0, 160000);

   TH2F *histSiPIDTraceIntegral = new TH2F("histSiPIDTraceIntegral", "histSiPIDTraceIntegral", 4000, 0, 90000, 4000, 0, 140000);
   TH2F *histSiPIDADCMax = new TH2F("histSiPIDADCMax", "histSiPIDADCMax", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax3He = new TH2F("histSiPIDADCMax3He", "histSiPIDADCMax3He", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax3HeGaggBe = new TH2F("histSiPIDADCMax3HeGaggBe", "histSiPIDADCMax3HeGaggBe", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax2H = new TH2F("histSiPIDADCMax2H", "histSiPIDADCMax2H", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMax1H = new TH2F("histSiPIDADCMax1H", "histSiPIDADCMax1H", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMaxGaggBe = new TH2F("histSiPIDADCMaxGaggBe", "histSiPIDADCMaxGaggBe", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMaxGaggLi = new TH2F("histSiPIDADCMaxGaggLi", "histSiPIDADCMaxGaggLi", 1000, 0, 4000, 1000, 0, 4000);
   TH1F *histSi1Cal = new TH1F("histSi1Cal", "histSi1Cal;ADCMaxFront;Counts", 1000, 0, 10000);
   TH1F *histSi2Cal = new TH1F("histSi2Cal", "histSi2Cal;ADCMaxFront;Counts", 1000, 0, 10000);

   TH1F *histSi1FrontCal = new TH1F("histSi1FrontCal", "histSi1FrontCal;E1_Front;Counts", 5000, 0, 50000);
   TH1F *histSi1BackCal = new TH1F("histSi1BackCal", "histSi1BackCal;E1_Back;Counts", 5000, 0, 50000);
   TH1F *histSi2FrontCal = new TH1F("histSi2FrontCal", "histSi2FrontCal;E2_Front;Counts", 5000, 0, 50000);
   TH1F *histSi2BackCal = new TH1F("histSi2BackCal", "histSi2BackCal;E2_Back;Counts", 5000, 0, 50000);

   TH1F *histSiMultiplicityFront1 = new TH1F("histSiMultiplicityFront1", "histSiMultiplicityFront1", 5, 0, 5);
   TH1F *histSiMultiplicityFront2 = new TH1F("histSiMultiplicityFront2", "histSiMultiplicityFront2", 5, 0, 5);

   TH1F *histGAGGCal = new TH1F("histGAGGCal", "histGAGGCal;E1;Counts", 5000, 0, 10000);
   TH1F *histGAGG1Cal = new TH1F("histGAGG1Cal", "histGAGG1Cal;ADCMaxFront;Counts", 5000, 0, 10000);
   TH1F *histGAGG2Cal = new TH1F("histGAGG2Cal", "histGAGG2Cal;ADCMaxFront;Counts", 5000, 0, 10000);
   TH1F *histGaggMultiplicity1 = new TH1F("histGaggMultiplicity1", "histGaggMultiplicity1", 25, 0, 25);
   TH1F *histGaggMultiplicity2 = new TH1F("histGaggMultiplicity2", "histGaggMultiplicity2", 16, 0, 16);
   TH2F *histGaggPIDADCMax = new TH2F("histGaggPIDADCMax", "histGaggPIDADCMax", 5000, 0, 10000, 4000, 0, 4000);
   TH2F *histGaggPIDADCMax3HeSiB = new TH2F("histGaggPIDADCMax3HeSiB", "histGaggPIDADCMax3HeSiB", 5000, 0, 10000, 4000, 0, 4000);
   TH2F *histGaggPIDADCMaxSiB = new TH2F("histGaggPIDADCMaxSiB", "histGaggPIDADCMaxSiB", 5000, 0, 10000, 4000, 0, 4000);

   // File input. There has to be a more compressed way to do this. As in, specify 3000 and then number of files
    std::vector runNums = {5139, 5140, 5141, 5142, 5143, 5144, 5145, 5146, 5147, 5148, 5149, 5150, 5151, 5152, 5153, 5154, 5155, 5156, 5157, 5158, 5159, 5160}; //Full calibration runs
    //std::vector runNums = {5139}; //Debug

   for (int runNum: runNums) {
      // Open the digitalization file and get the TTree.
      TString unpackFileName = TString::Format("/home/astinson/e535rawdata/UnpackerTestOutput/run_%04d.root", runNum);
      TFile *unpackFile = new TFile(unpackFileName, "READ");
      TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
      int nUnpackEvents = unpackTree->GetEntries();
      std::cout << " Number of unpacked events in run " << runNum << ": " << nUnpackEvents << std::endl;
      int nEventsWith2Tracks = 0;
      int nEventsWith3He = 0;
      // Creare the TTreeReader to read the AtTrackingEvents and simulation.
      TTreeReader unpackReader("cbmsim", unpackFile);
      TTreeReaderValue<TClonesArray> siArray(unpackReader, "AtSiEvent");
      TTreeReaderValue<TClonesArray> gaggArray(unpackReader, "AtGaggEvent");
      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");

      // Loop over events.
      for (int i = 0; i < nUnpackEvents; i++) {
      //for (int i = 0; i < nGaggEvents; i++) {
         unpackReader.Next();

         // Check the Si data first.
         AtSiEvent *siEvent = (AtSiEvent *)siArray->At(0);

         Int_t multiplicityFront1 = siEvent->GetMultiplicityFront1(); //Hits detected in each layer (1vs2)
         Int_t multiplicityFront2 = siEvent->GetMultiplicityFront2();

         histSiMultiplicityFront1->Fill(multiplicityFront1);
         histSiMultiplicityFront2->Fill(multiplicityFront2);
	 
	 Double_t maxADCFront1 = siEvent->GetADCMaxFront1(0);
	 Double_t maxADCFront2 = siEvent->GetADCMaxFront2(0);

	 histSi1Cal->Fill(maxADCFront1); //Plotting raw ADC to 1D hist
	 histSi2Cal->Fill(maxADCFront2);

         //if (multiplicityFront1 != 1 || multiplicityFront2 != 1) continue; //Only single hits allowed

         Double_t EFront1 = siEvent->GetEFront1(0); //CHECK energy deposited in the Si
         Double_t EFront2 = siEvent->GetEFront2(0);
	 Double_t EBack1 = siEvent->GetEBack1(0);
	 Double_t EBack2 = siEvent->GetEBack2(0);

	 histSi1FrontCal->Fill(EFront1);
	 histSi1BackCal->Fill(EBack1);
	 histSi2FrontCal->Fill(EFront2);
	 histSi2BackCal->Fill(EBack2);
         histSiPIDTraceIntegral->Fill(EFront2, EFront1);
         histSiPIDADCMax->Fill(maxADCFront2, maxADCFront1);
	 //histSi1Cal->Fill(maxADCFront1); //Plotting raw ADC to 1D hist
	 //histSi2Cal->Fill(maxADCFront2);

         // Check the GAGG data
         AtGaggEvent *gaggEvent = (AtGaggEvent *)gaggArray->At(0);

         Int_t multiGagg1 = gaggEvent->GetMultiplicity1();
         Int_t multiGagg2 = gaggEvent->GetMultiplicity2();

         histGaggMultiplicity1->Fill(multiGagg1);
         histGaggMultiplicity2->Fill(multiGagg2);

	 Double_t EG1, EG2;
	 EG1 = gaggEvent->GetE1(0);

	 histGAGGCal->Fill(EG1);
	 /*Double_t ADCG1, ADCG2;
	 ADCG1 = gaggEvent->GetADCMax1(0);
	 ADCG2 = gaggEvent->GetADCMax2(0);

	 histGAGG1Cal->Fill(ADCG1);
	 histGAGG2Cal->Fill(ADCG2);*/

         //std::cout << " Multiplicity GAGG1: " << multiGagg1 <<std::endl;

         //if (multiGagg1 != 1) continue;
         Double_t maxADCGagg1, maxADCGagg2;
         Double_t gaggADC1{-1}, gaggADC2{-1};
         Double_t totalMaxADCGagg{};
         for(int imul = 0; imul < multiGagg1; imul++) {
            maxADCGagg1 = gaggEvent->GetADCMax1(imul);
            gaggADC1 += maxADCGagg1;
            totalMaxADCGagg += maxADCGagg1;
            //std::cout << " ADC entry = " << maxADCGagg1 << std::endl;
         }
	      for(int imul = 0; imul < multiGagg2; imul++) {
            maxADCGagg2 = gaggEvent->GetADCMax2(imul); 
            gaggADC2 += maxADCGagg2;
            totalMaxADCGagg += maxADCGagg2;
         }

	 histGAGGCal->Fill(totalMaxADCGagg); //Here Amber for GAGG
	 histGaggPIDADCMax->Fill(totalMaxADCGagg, maxADCFront2);


         // First, we obtain some rough kinematics just by using the AtPatternEvent.
         AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
         if (!patternEvent) continue;

         // We want to focus on events with 2 or less tracks for now.
         auto &tracks = patternEvent->GetTrackCand();
         //if (tracks.size() > 4) continue;
         if (tracks.size() > 2) continue;
         //if (tracks.size() > 1) continue;
	 if (tracks.size()==2) nEventsWith2Tracks++;
	 int trackIndex = 0;
	 double thetaLABArray[2] = {0.0,0.0};
	 double estimatedKineEArray[2] = {0.0,0.0};
	 bool isHe3Array[2] = {false,false};
         // Iterate over AtTracks and extract their kinematics.
         for (auto &track: tracks) {
            bool isPunchThrough = punchThroughChecker.IsPunchThrough(&track);
            if (isPunchThrough) continue;

            auto *pattern = track.GetPattern();

            auto firstPoint = track.GetFirstPoint();
            auto lastPoint = track.GetLastPoint();
	    //cut penetrate particles
	    //if (lastPoint.Z()>900) continue;
	    //std::cout<<lastPoint.Z()<<std::endl;
            double roughRangeEstimation = pattern->DistanceAlongPattern(lastPoint, firstPoint);

            double trackThetaLAB = track.GetGeoTheta() * 180 / TMath::Pi();
            double trackPhi = track.GetGeoPhi() * 180 / TMath::Pi();

            double smallPadCharge{};
            double bigPadCharge{};
            auto braggCurvePairs = track.GetBraggCurveValues();
            auto &hits = track.GetHitArray();
            double rangeInSmallPads{};
	         for (auto &hit: hits) {
               int padNum = hit->GetPadNum();
               int sizeID = map->GetPadSize(padNum);
               if (sizeID == 1) {
                  bigPadCharge += hit->GetCharge();
                  continue;
               }
               smallPadCharge += hit->GetCharge();
               double currentRangeInSmallPads = pattern->DistanceAlongPattern(hit->GetPosition(), firstPoint);
               if (currentRangeInSmallPads > rangeInSmallPads)
                  rangeInSmallPads = currentRangeInSmallPads;
            }
            double dEdx = smallPadCharge / rangeInSmallPads;

            double rangeInBigPads = roughRangeEstimation - rangeInSmallPads;
            bool reachedBigPads = true;
            if (rangeInBigPads / roughRangeEstimation < 0.05)
               reachedBigPads = false;

            double estimatedKinE{0.1};
            if (cutPID3He->IsInside(roughRangeEstimation, dEdx)) {
               while (eLossModelC3D8_3He->GetRange(estimatedKinE) < roughRangeEstimation)
                  estimatedKinE += 0.01;
	    } else if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx)) {
               while (eLossModelC3D8_d->GetRange(estimatedKinE) < roughRangeEstimation)
                  estimatedKinE += 0.01;
	    } else if (cutPIDproton->IsInside(roughRangeEstimation, dEdx)) {
               while (eLossModelC3D8_p->GetRange(estimatedKinE) < roughRangeEstimation)
                  estimatedKinE += 0.01;
            } else {
	      estimatedKinE = 0;
            }

            histRangeVThetaLAB->Fill(trackThetaLAB, roughRangeEstimation);
            histdEdxVTotalRange->Fill(roughRangeEstimation, dEdx);
            histEstimatedKinEVThetaLABTotal->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx))
               histTrackThetaLABRange3He->Fill(roughRangeEstimation, trackThetaLAB);

            if ((cutPID3He->IsInside(roughRangeEstimation, dEdx))
		&& cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2)
		&& cutSiBe->IsInside(maxADCFront2, maxADCFront1)){
	      histTrackThetaLABRange3HeSiGaggBe->Fill(roughRangeEstimation, trackThetaLAB);};

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB3He->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx)
		&& cutSiBe->IsInside(maxADCFront2, maxADCFront1)){
	      histEstimatedKinEVThetaLAB3HeSiBe->Fill(trackThetaLAB, estimatedKinE);};

            if (cutSiBe->IsInside(maxADCFront2, maxADCFront1))
	      histEstimatedKinEVThetaLABSiBe->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx)
		&& cutSiBe->IsInside(maxADCFront2, maxADCFront1)
		&& cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2)){
	      histEstimatedKinEVThetaLAB3HeSiGaggBe->Fill(trackThetaLAB, estimatedKinE);
	      if(tracks.size()==1){
		histEstimatedKinEVThetaLAB3HeSiGaggBetrack1->Fill(trackThetaLAB, estimatedKinE);
	      }
	      if(tracks.size()==2){
		histEstimatedKinEVThetaLAB3HeSiGaggBetrack2->Fill(trackThetaLAB, estimatedKinE);
	      }
	    };

            if (cutSiBe->IsInside(maxADCFront2, maxADCFront1)
		&& cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2)){
	      histEstimatedKinEVThetaLABSiGaggBe->Fill(trackThetaLAB, estimatedKinE);};

            if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB2H->Fill(trackThetaLAB, estimatedKinE);

            if (cutPIDproton->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB1H->Fill(trackThetaLAB, estimatedKinE);
	    
	    if (cutSiB->IsInside(maxADCFront2, maxADCFront1))
	      histEstimatedKinEVThetaLABSiB->Fill(trackThetaLAB, estimatedKinE);

	    if (cutSi13B->IsInside(maxADCFront2, maxADCFront1))
	      histEstimatedKinEVThetaLABSi13B->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx))
	      histSiPIDADCMax3He->Fill(maxADCFront2, maxADCFront1);

            if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx))
	      histSiPIDADCMax2H->Fill(maxADCFront2, maxADCFront1);

	    if (cutPIDproton->IsInside(roughRangeEstimation, dEdx))
	      histSiPIDADCMax1H->Fill(maxADCFront2, maxADCFront1);

	    if (cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2))
	      histSiPIDADCMaxGaggBe->Fill(maxADCFront2, maxADCFront1);

            if ((cutPID3He->IsInside(roughRangeEstimation, dEdx))
		&& (cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2))){
	      histSiPIDADCMax3HeGaggBe->Fill(maxADCFront2, maxADCFront1);};

	    if (cutGaggLi->IsInside(totalMaxADCGagg, maxADCFront2))
	      histSiPIDADCMaxGaggLi->Fill(maxADCFront2, maxADCFront1);

	    if (cutSiC->IsInside(maxADCFront2, maxADCFront1))
	      histEstimatedKinEVThetaLABSiC->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx)
		&& cutSiB->IsInside(maxADCFront2, maxADCFront1)){
	      histGaggPIDADCMax3HeSiB->Fill(totalMaxADCGagg, maxADCFront2);};

	    if (cutSiB->IsInside(maxADCFront2, maxADCFront1))
	      histGaggPIDADCMaxSiB->Fill(totalMaxADCGagg, maxADCFront2);

            if (!reachedBigPads)
               histESmallVTotalRange->Fill(roughRangeEstimation, smallPadCharge);
            else
               histEBigVBigRange->Fill(rangeInBigPads, bigPadCharge);

	    if(tracks.size()==2 && trackIndex < 2){
	      thetaLABArray[trackIndex] = trackThetaLAB;
	      estimatedKineEArray[trackIndex] = estimatedKinE;
	      isHe3Array[trackIndex] = cutPID3He->IsInside(roughRangeEstimation, dEdx);
	      trackIndex++;
	    }
         }
         if(tracks.size()==2 && trackIndex==2){
	   histThetaLABThetaLAB->Fill(thetaLABArray[0], thetaLABArray[1]);
	   if(thetaLABArray[0]<=thetaLABArray[1]){
	     if(thetaLABArray[0]<15 && thetaLABArray[1]<50){
	       histThetaLABThetaLABKine->Fill(thetaLABArray[1], thetaLABArray[0]);
	       histEstimatedKinEVThetaLABkinett->Fill(thetaLABArray[0],estimatedKineEArray[0]);
	     }
	   }
	   else if(thetaLABArray[0]>thetaLABArray[1]){
	     if(thetaLABArray[0]<50 && thetaLABArray[1]<15){
	       histThetaLABThetaLABKine->Fill(thetaLABArray[0], thetaLABArray[1]);
	       histEstimatedKinEVThetaLABkinett->Fill(thetaLABArray[1],estimatedKineEArray[1]);
	     }
	   }
	   bool passSiBe = cutSiBe->IsInside(maxADCFront2, maxADCFront1);
	   bool passGaggBe = cutGaggBe->IsInside(totalMaxADCGagg, maxADCFront2);
	   if(passSiBe && passGaggBe){
	     if(isHe3Array[0] && !isHe3Array[1]){
               histThetaLABThetaLAB3HeSiGaggBe->Fill(thetaLABArray[0], thetaLABArray[1]);
	       nEventsWith3He++ ;
	     }
	     else if(!isHe3Array[0] && isHe3Array[1]){
               histThetaLABThetaLAB3HeSiGaggBe->Fill(thetaLABArray[1], thetaLABArray[0]);
	       nEventsWith3He++;
	     }
	   }
         }
      }
      //      std::cout << "Number of 2 tracks events in run" << runNum << ":" << nEventsWith2Tracks << std::endl;
      //      std::cout << "Number of events with 3He in run" << runNum << ":" << nEventsWith3He << std::endl;
      // Close files.
      unpackFile->Close();
      //gaggFile->Close();
   }

   outfile->Write();

   // Kinematic lines.
   TGraph *kine_dd = ReadKinematics("./kineFiles/17N_dd_gs.txt");
   TGraph *kine_d3He = ReadKinematics("./kineFiles/13B_d3He_gs.txt");
   TGraph *kine_d3HeEx2_2 = ReadKinematics("./kineFiles/13B_d3He_Ex2_2.txt");
   TGraph *kine_d3HeEx2_7 = ReadKinematics("./kineFiles/13B_d3He_Ex2_7.txt");
   TGraph *kine_d3He_tt = ReadKinematics("./kineFiles/13B_d3He_tt.txt");
   TGraph *kine_12C12C = ReadKinematics("./kineFiles/17N_12C12C_gs.txt");

   
   // Draw histograms in TCanvas.
   outfile->Close();
}

TGraph* ReadKinematics(TString kineFile)
{
   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()){
      while (!kineStr->eof()){
	//         *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >>
	//                     ThetaLabSca[numKin] >> EnerLabSca[numKin];
         *kineStr >> ThetaLabRec[numKin] >> EnerLabRec[numKin];
         numKin++;
      }
   } else if (kineStr->fail())
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   return kine;
}

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{

   // in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables

   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) -
                 (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) /
                   (2 * pow(m2, 2)) +
                s + u - pow(m2, 2));
   Ex = m4_ex - m4;

   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;

   // for inverse kinematics Note: this angle corresponds to the recoil
    theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   /*theta_cm = acos((pow(s, 2) + s * (2 * u - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m4_ex, 2) - pow(m3, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m4_ex, 2), pow(m3, 2))));*/

   theta_cm = theta_cm * TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}
