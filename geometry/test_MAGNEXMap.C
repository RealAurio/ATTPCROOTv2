void test_MAGNEXMap()
{

   TString mapFile = "/home/aurio/research/NUMEN/NUMEN_data/channel2pad_2.txt";
   TString parFile = "/home/aurio/research/NUMEN/NUMEN_data/testParFile.txt";


   gSystem->Load("libAtTpcMap.so");
   AtMAGNEXMap *magnex = new AtMAGNEXMap(parFile, mapFile);
   magnex->GeneratePadPlane();

   TH2Poly *padplane = magnex->GetPadPlane();

   for (auto i = 0; i < 20; ++i) {

      Int_t bin = padplane->Fill(i + 0.5, -i - 0.5, i * 100);
      std::cout << " Bin : " << bin << "\n";
   }

   padplane->Draw();
   //gPad->Update();
}
