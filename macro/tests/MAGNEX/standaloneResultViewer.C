void standaloneResultViewer(int runNumber = 38)
{
   // Define the histograms to be filled.
   TH1F *xRansacEntranceHist = new TH1F("xRansacEntranceHist", "xRansacEntranceHist", 600, 0, 300);
   TH1F *yRansacEntranceHist = new TH1F("yRansacEntranceHist", "yRansacEntranceHist", 800, -100, 300);
   TH1F *zRansacEntranceHist = new TH1F("zRansacEntranceHist", "zRansacEntranceHist", 120, -10, 110);

   TH1F *thetaHist = new TH1F("thetaHist", "thetaHist", 2880, -180, 180);
   TH1F *phiHist = new TH1F("phiHist", "phiHist", 2280, -180, 180);

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

   // Open the results file and load the data.
   Int_t eventID{}, trackID{};
   Long64_t firstTrackTimeStamp{};
   Double_t xRansacEntrance{}, yRansacEntrance{}, zRansacEntrance{}, theta{}, phi{};
   Double_t xRansacRow[5], yRansacRow[5], zRansacRow[5];

   TFile *resultFile = new TFile(TString::Format("/home/aurio/research/NUMEN/ATTPCROOTv2_stuff/run%03d.root", runNumber), "READ");
   TTree *anatree = (TTree *)resultFile->Get("anatree");
   Int_t nEntries = anatree->GetEntries();
   std::cout << " Number of entries : " << nEntries << std::endl;

   anatree->SetBranchAddress("eventID", &eventID);
   anatree->SetBranchAddress("trackID", &trackID);
   anatree->SetBranchAddress("firstTrackTimeStamp", &firstTrackTimeStamp);
   anatree->SetBranchAddress("xRansacEntrance", &xRansacEntrance);
   anatree->SetBranchAddress("yRansacEntrance", &yRansacEntrance);
   anatree->SetBranchAddress("zRansacEntrance", &zRansacEntrance);
   anatree->SetBranchAddress("theta", &theta);
   anatree->SetBranchAddress("phi", &phi);
   anatree->SetBranchAddress("xRansacRow", xRansacRow);
   anatree->SetBranchAddress("yRansacRow", yRansacRow);
   anatree->SetBranchAddress("zRansacRow", zRansacRow);

   // File where to write the TS, angle, event and track ID.
   std::ofstream outputFile;
   outputFile.open(TString::Format("ATTPCROOTv2_RANSACangles_run%03d.txt", runNumber).Data());

   // Iterate over entries.
   std::cout << "Entry:" << std::endl;
   for (Int_t i = 0; i < nEntries; i++) {
      if (i % 10000 == 0)
         std::cout << "\x1b[1A" << "\x1b[2K" << "Entry: " << i << std::endl;

      anatree->GetEntry(i);

      xRansacEntranceHist->Fill(xRansacEntrance);
      yRansacEntranceHist->Fill(yRansacEntrance);
      zRansacEntranceHist->Fill(zRansacEntrance);

      thetaHist->Fill(theta);
      phiHist->Fill(phi);

      for (int iRow = 0; iRow < 5; iRow++) {
         xRansacRowHist[iRow]->Fill(xRansacRow[iRow]);
         yRansacRowHist[iRow]->Fill(yRansacRow[iRow]);
         zRansacRowHist[iRow]->Fill(zRansacRow[iRow]);
      }

      outputFile << " Event: " << eventID << " | Track: " << trackID << "\n";
      outputFile << "   firstTrackTimeStamp = " << firstTrackTimeStamp << " ps\n";
      outputFile << "   theta = " << theta << " deg\n";
      outputFile << "   phi = " << phi << " deg\n";
      outputFile << std::endl;

   }
   // Close the files.
   resultFile->Close();
   outputFile.close();

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
}
