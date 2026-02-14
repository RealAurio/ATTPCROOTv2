void histogrammer(int runNumber = 38)
{
   using XYPoint = ROOT::Math::XYPoint;
   using XYZPoint = ROOT::Math::XYZPoint;
   using XYZVector = ROOT::Math::XYZVector;

   // MAGNEX pad map.
   AtMAGNEXMap *map = new AtMAGNEXMap("/home/aurio/research/NUMEN/NUMEN_data/testParFile.txt");
   map->GeneratePadPlane();

   // Define the histograms to be filled.
   TH1F *xRansacEntranceHist = new TH1F("xRansacEntranceHist", "xRansacEntranceHist", 600, 0, 300);
   TH1F *yRansacEntranceHist = new TH1F("yRansacEntranceHist", "yRansacEntranceHist", 800, -100, 300);
   TH1F *zRansacEntranceHist = new TH1F("zRansacEntranceHist", "zRansacEntranceHist", 120, -10, 110);

   TH1F *thetaHist = new TH1F("thetaHist", "thetaHist", 2880, -180, 180);
   TH1F *phiHist = new TH1F("phiHist", "phiHist", 2280, -180, 180);

   TH1F *multiplicityHist = new TH1F("multiplicityHist", "multiplicityHist", 10, -0.5, 9.5);

   std::vector<TH1F *> xRansacRowHist = {new TH1F("xRansacRow0", "xRansacRow0", 600, 0, 300),
                                         new TH1F("xRansacRow1", "xRansacRow1", 600, 0, 300),
                                         new TH1F("xRansacRow2", "xRansacRow2", 600, 0, 300),
                                         new TH1F("xRansacRow3", "xRansacRow3", 600, 0, 300),
                                         new TH1F("xRansacRow4", "xRansacRow4", 600, 0, 300)};

   std::vector<TH1F *> yRansacRowHist = {new TH1F("yRansacRow0", "yRansacRow0", 800, -100, 300),
                                         new TH1F("yRansacRow1", "yRansacRow1", 800, -100, 300),
                                         new TH1F("yRansacRow2", "yRansacRow2", 800, -100, 300),
                                         new TH1F("yRansacRow3", "yRansacRow3", 800, -100, 300),
                                         new TH1F("yRansacRow4", "yRansacRow4", 800, -100, 300)};

   std::vector<TH1F *> zRansacRowHist = {new TH1F("zRansacRow0", "zRansacRow0", 120, -10, 110),
                                         new TH1F("zRansacRow1", "zRansacRow1", 120, -10, 110),
                                         new TH1F("zRansacRow2", "zRansacRow2", 120, -10, 110),
                                         new TH1F("zRansacRow3", "zRansacRow3", 120, -10, 110),
                                         new TH1F("zRansacRow4", "zRansacRow4", 120, -10, 110)};

   // Related to the calculation of the source point. It probably does not apply to every case.
   TH1F *xSourceHist = new TH1F("xSourceHist", "xSourceHist", 1100, -100, 1000);
   TH1F *zSourceHist = new TH1F("zSourceHist", "zSourceHist", 100, -200, 200);
   TH2F *sourceHist = new TH2F("sourceHist", "sourceHist", 250, 0, 500, 250, -300, 200);
   Double_t x0Prev{}, z0Prev{};
   Double_t xPrimePrev{}, zPrimePrev{};
   Bool_t firstTrack{true};

   // Create a regular root file that can be opened with stand-alone root.
   Int_t eventID{}, trackID{};
   Long64_t firstTrackTimeStamp{};
   Double_t xRansacEntrance{}, yRansacEntrance{}, zRansacEntrance{}, theta{}, phi{};
   Double_t xRansacRow[5], yRansacRow[5], zRansacRow[5];

   TFile *outFile = new TFile(TString::Format("/home/aurio/research/NUMEN/ATTPCROOTv2_stuff/run%03d.root", runNumber), "RECREATE");
   TTree *anatree = new TTree("anatree", "TTree with values reconstructed by ATTPCROOTv2");
   anatree->Branch("eventID", &eventID);
   anatree->Branch("trackID", &trackID);
   anatree->Branch("firstTrackTimeStamp", &firstTrackTimeStamp);
   anatree->Branch("xRansacEntrance", &xRansacEntrance);
   anatree->Branch("yRansacEntrance", &yRansacEntrance);
   anatree->Branch("zRansacEntrance", &zRansacEntrance);
   anatree->Branch("theta", &theta);
   anatree->Branch("phi", &phi);
   anatree->Branch("xRansacRow", xRansacRow, "xRansacRow[5]/D");
   anatree->Branch("yRansacRow", yRansacRow, "yRansacRow[5]/D");
   anatree->Branch("zRansacRow", zRansacRow, "zRansacRow[5]/D");

   // Initialize the FairRun, and load the output of the recontruct.C macro.
   FairRunAna *run = new FairRunAna();

   TFile *file = new TFile(TString::Format("/home/aurio/research/NUMEN/NUMEN_output/dataBrazil/reconstructed_%03d.root", runNumber), "READ");
   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtHitClusterEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");

   // Iterate over events.
   std::cout << "Event:" << std::endl;
   //for (Int_t i = 0; i < nEvents; i++) {
   for (Int_t i = 0; i < 200000; i++) {
      if (i % 10000 == 0)
         std::cout << "\x1b[1A" << "\x1b[2K" << "Event: " << i << std::endl;

      Reader.Next();
      AtEvent *event = (AtEvent *)eventHArray->At(0);
      AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

      // If for some reason there is not a valid event, skip.
      if (!event || !patternEvent)
         continue;

      // Access the tracks in each event.
      auto &tracks = patternEvent->GetTrackCand();
      multiplicityHist->Fill(tracks.size());

      //if (tracks.size() > 2)
         //std::cout << "Event " << i << " has " << tracks.size() << " tracks." << std::endl;

      eventID = i;
      trackID = -1;
      for (auto &track : tracks) {
         trackID++;

         firstTrackTimeStamp = track.GetHitArray()[0]->GetTimeStamp();

         // Access the fitted pattern line in each track.
         const AtPatterns::AtPatternLine *patternLine = dynamic_cast<const AtPatterns::AtPatternLine *>(track.GetPattern());

         // Position in the entrance plane of the TPC. The z component should be 0 by construction.
         XYZPoint entrancePoint = patternLine->GetPointAt(0);
         xRansacEntrance = entrancePoint.X();
         yRansacEntrance = entrancePoint.Y();
         zRansacEntrance = entrancePoint.Z();

         // Direction vector of the track. The z component should be 1 or -1 by construction, and the x and
         // y components the tangents of the entrance angles.
         XYZVector directionVector = patternLine->GetDirection();
         theta = TMath::ATan2(directionVector.X(), directionVector.Z()) * 180. / TMath::Pi();
         phi = TMath::ATan2(directionVector.Y(), directionVector.Z()) * 180. / TMath::Pi();

         // Fill the histograms.
         //if (y < 35 || 65 < y)
            //continue;

         xRansacEntranceHist->Fill(xRansacEntrance);
         yRansacEntranceHist->Fill(yRansacEntrance);
         zRansacEntranceHist->Fill(zRansacEntrance);

         thetaHist->Fill(theta);
         phiHist->Fill(phi);

         for (int iRow = 0; iRow < 5; iRow++) {
            XYPoint point = map->CalcPadCenter(map->PadID(0, iRow));
            Double_t zRow = point.Y();

            XYZPoint rowPoint = patternLine->GetPointAt(zRow);
            xRansacRow[iRow] = rowPoint.X();
            yRansacRow[iRow] = rowPoint.Y();
            zRansacRow[iRow] = rowPoint.Z();

            xRansacRowHist[iRow]->Fill(xRansacRow[iRow]);
            yRansacRowHist[iRow]->Fill(yRansacRow[iRow]);
            zRansacRowHist[iRow]->Fill(zRansacRow[iRow]);
         }

         // Related to the calculation of the source point. It probably does not apply to every case.
         if (!firstTrack) {
            Double_t zSource = (entrancePoint.X() - x0Prev + xPrimePrev / zPrimePrev * z0Prev - directionVector.X() / directionVector.Z() * entrancePoint.Z()) / (xPrimePrev / zPrimePrev - directionVector.X() / directionVector.Z());
            Double_t xSource = x0Prev + xPrimePrev * (zSource - z0Prev) / zPrimePrev;

            xSourceHist->Fill(xSource);
            zSourceHist->Fill(zSource);
            sourceHist->Fill(xSource, zSource);
         }
         firstTrack = false;

         x0Prev = entrancePoint.X();
         z0Prev = entrancePoint.Z();
         xPrimePrev = directionVector.X();
         zPrimePrev = directionVector.Z();

         // Fill TTree.
         if (trackID > -1)
            anatree->Fill();
      }
   }
   // Write TTree in the output file.
   outFile->cd();
   anatree->Write();

   // Close the files.
   file->Close();
   outFile->Close();

   // Draw the histograms.
   TCanvas *cXRansacEntrance = new TCanvas();
   xRansacEntranceHist->Draw();
   xRansacEntranceHist->GetXaxis()->SetTitle("x [mm]");

   TCanvas *cYRansacEntrance = new TCanvas();
   yRansacEntranceHist->Draw();
   yRansacEntranceHist->GetXaxis()->SetTitle("y [mm]");

   TCanvas *cZRansacEntrance = new TCanvas();
   zRansacEntranceHist->Draw();
   zRansacEntranceHist->GetXaxis()->SetTitle("z [mm]");

   TCanvas *cTheta = new TCanvas();
   thetaHist->Draw();
   thetaHist->GetXaxis()->SetTitle("#theta [deg]");

   TCanvas *cPhi = new TCanvas();
   phiHist->Draw();
   phiHist->GetXaxis()->SetTitle("#phi [deg]");

   TCanvas *cRansacPerRow = new TCanvas();
   cRansacPerRow->Divide(5, 2);
   std::vector<TCanvas *> cZRansacRow = {};
   for (int iRow = 0; iRow < 5; iRow++) {
      cRansacPerRow->cd(iRow + 1);
      xRansacRowHist[iRow]->Draw();
      xRansacRowHist[iRow]->GetXaxis()->SetTitle("x [mm]");

      cRansacPerRow->cd(iRow + 6);
      yRansacRowHist[iRow]->Draw();
      yRansacRowHist[iRow]->GetXaxis()->SetTitle("y [mm]");

      cZRansacRow.push_back(new TCanvas());
      zRansacRowHist[iRow]->Draw();
      zRansacRowHist[iRow]->GetXaxis()->SetTitle("z [mm]");
   }

   TCanvas *cMultiplicity = new TCanvas();
   multiplicityHist->Draw();
   multiplicityHist->GetXaxis()->SetTitle("multiplicity");

   TCanvas *cXSource = new TCanvas();
   xSourceHist->Draw();
   xSourceHist->GetXaxis()->SetTitle("x [mm]");

   TCanvas *cZSource = new TCanvas();
   zSourceHist->Draw();
   zSourceHist->GetXaxis()->SetTitle("z [mm]");

   TCanvas *cSource = new TCanvas();
   sourceHist->Draw("zcol");
   sourceHist->GetXaxis()->SetTitle("x [mm]");
   sourceHist->GetXaxis()->SetTitle("z [mm]");
}
