#include "AtBeamIonGenerator.h"

#include "AtVertexPropagator.h"

#include <FairIon.h>
#include <FairPrimaryGenerator.h>
#include <FairRunSim.h>

#include <TRandom.h>

#include <cmath>

constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";
constexpr auto cBLUE = "\033[1;34m";

AtBeamIonGenerator::AtBeamIonGenerator(TString ionName, int pdgCode, int Z, int A, int chargeState, double mass,
                                       double Ex, double kineticEnergy, double sigmaKE)
   : fIonName(ionName), fZ(Z), fA(A), fChargeState(chargeState), fMass(mass), fEx(Ex)
{
   if (IonAlreadyExists(ionName)) {
      LOG(warning) << ionName.Data() << " ion already exists! Using the already defined one.";
   } else {
      LOG(info) << " Adding new ion named " << ionName.Data() << " to the simulation.";
      LOG(info) << "  Z = " << Z << " | A = " << A << " | chargeState = " << chargeState;
      LOG(info) << "  mass = " << mass << " MeV | Ex = " << Ex << " MeV";
      fIon = new FairIon(ionName.Data(), Z, A, chargeState, Ex / 1000, mass / 1000);

      if (!AddIonToFairRunSim())
         return;
   }

   fPDGCode = pdgCode;

   ConfigureVertexPropagator();

   if (kineticEnergy < 0)
      return;
   if (sigmaKE < 0)
      SetBeamEnergy(kineticEnergy);
   else
      SetBeamEnergy(kineticEnergy, sigmaKE);
}

AtBeamIonGenerator::AtBeamIonGenerator(TString ionName, int pdgCode, double kineticEnergy, double sigmaKE)
   : fIonName(ionName)
{
   if (!IonAlreadyExists(ionName)) {
      LOG(fatal) << ionName.Data() << " ion was not found! Please use the other initializer instead!";
      return;
   }
   LOG(info) << "Using the already defined " << ionName.Data() << " ion as requested.";

   fZ = fIon->GetZ();
   fA = fIon->GetA();
   fChargeState = fIon->GetQ();
   fMass = fIon->GetMass() * 1000;
   fEx = fIon->GetExcEnergy() * 1000;
   fPDGCode = pdgCode;

   ConfigureVertexPropagator();

   if (kineticEnergy < 0)
      return;
   if (sigmaKE < 0)
      SetBeamEnergy(kineticEnergy);
   else
      SetBeamEnergy(kineticEnergy, sigmaKE);
}

AtBeamIonGenerator::~AtBeamIonGenerator()
{
   if (fFileDistribution)
      fFileDistribution->Close();
}

bool AtBeamIonGenerator::IonAlreadyExists(TString ionName)
{
   // Extract the ions from the FairRunSim.
   FairRunSim *run = FairRunSim::Instance();
   if (!run) {
      LOG(error) << "No FairRunSim instantised! Can not get ion list from it!";
      return false;
   }
   TObjArray *UserIons = run->GetUserDefIons();

   // Check if there is an ion with the name.
   fIon = dynamic_cast<FairIon *>(UserIons->FindObject(ionName));
   if (fIon)
      return true;
   return false;
}

bool AtBeamIonGenerator::AddIonToFairRunSim()
{
   FairRunSim *run = FairRunSim::Instance();
   if (!run) {
      LOG(error) << "No FairRunSim instantised! Can not add ion to it!";
      return false;
   }
   run->AddNewIon(fIon);
   return true;
}

void AtBeamIonGenerator::ConfigureVertexPropagator()
{
   // TO-DO: ADD MORE CONFIGURATION OPTIONS.
   AtVertexPropagator::Instance()->SetBeamMass(fIon->GetMass() + fIon->GetExcEnergy());
}

void AtBeamIonGenerator::SetBeamSource(double z, double x, double y)
{
   fSourceX = x;
   fSourceY = y;
   fSourceZ = z;
   fBeamProfileType = DistributionType::kDeltaValue;

   LOG(info) << " The source of the beam ion has been set to (" << fSourceX << ", " << fSourceY << ", " << fSourceZ
             << ") [mm].";
}

void AtBeamIonGenerator::SetBeamSource(double z, double x, double y, double beamRadius)
{
   SetBeamSource(z, x, y);

   if (beamRadius < 0)
      LOG(warning) << " You are setting negative beam radius for the beam profile! The simulation code should still "
                      "work fine, but maybe you did not intend to set that value...";
   fBeamRadius = beamRadius;
   fBeamProfileType = DistributionType::kUniform;

   LOG(info) << " Added a beam profile radius of " << fBeamRadius << " [mm].";
}

void AtBeamIonGenerator::SetBeamSourceGaussian(double z, double x, double y, double beamRadius)
{
   SetBeamSource(z, x, y, beamRadius);
   fBeamProfileType = DistributionType::kGaussian;

   LOG(info) << " Setting beam profile to gaussian mode instead of uniform mode.";
}

void AtBeamIonGenerator::SetBeamSource(double z, TString pathToDistributionFile)
{
   fSourceZ = z;

   if (!fFileDistribution)
      fFileDistribution = std::make_unique<TFile>(pathToDistributionFile.Data(), "READ");

   if (fFileDistribution->IsZombie()) {
      LOG(error) << " Tried to pass beam profile histogram, but the file " << pathToDistributionFile.Data()
                 << " could not be open!";
      fBeamProfileType = DistributionType::kUndefined;
      return;
   }

   fProfileDistribution = dynamic_cast<TH2F *>(fFileDistribution->Get("profileDistribution"));
   if (!fProfileDistribution) {
      LOG(error) << " Beam profile histogram named profileDistribution not found on the file "
                 << pathToDistributionFile.Data() << "!";
      fBeamProfileType = DistributionType::kUndefined;
      return;
   }
   fProfileDistribution->SetDirectory(nullptr);
   LOG(info) << " Beam profile histogram named profileDistribution found on the file " << pathToDistributionFile.Data()
             << ". Setting of beam profile succesful!";
   LOG(info) << " This beam profile histogram was set to correspond with a position z = " << fSourceZ
             << " [mm].";
   fBeamProfileType = DistributionType::kHistogram;
}

void AtBeamIonGenerator::SetBeamDirection(double thetaXZ, double thetaYZ)
{
   fThetaXZ = thetaXZ;
   fThetaYZ = thetaYZ;
   fBeamAnglesType = DistributionType::kDeltaValue;

   LOG(info) << " The beam direction has been set so that the non-dispersive angle (XZ) is " << fThetaXZ
             << " [deg], and the dispersive (YZ) is " << fThetaYZ << " [deg].";
}

void AtBeamIonGenerator::SetBeamDirection(double thetaXZ, double thetaYZ, double sigmaXZ, double sigmaYZ)
{
   SetBeamDirection(thetaXZ, thetaYZ);

   if (sigmaXZ < 0 || sigmaYZ < 0)
      LOG(warning) << " You are setting negative standard deviations for the angle distributions! The simulation code "
                      "should still work fine, but maybe you did not intend to set these values...";
   fSigmaXZ = sigmaXZ;
   fSigmaYZ = sigmaYZ;
   fBeamAnglesType = DistributionType::kGaussian;

   LOG(info) << " Added a standard deviation of the non-dispersive angle (XZ) of " << fSigmaXZ << " [deg], and "
             << fSigmaYZ << " [deg] for the dispersive one (YZ).";
}

void AtBeamIonGenerator::SetBeamDirection(TString pathToDistributionFile)
{
   if (!fFileDistribution)
      fFileDistribution = std::make_unique<TFile>(pathToDistributionFile.Data(), "READ");

   if (fFileDistribution->IsZombie()) {
      LOG(error) << " Tried to pass beam angular distribution histogram, but the file " << pathToDistributionFile.Data()
                 << " could not be opened!";
      fBeamAnglesType = DistributionType::kUndefined;
      return;
   }

   fAnglesDistribution = dynamic_cast<TH2F *>(fFileDistribution->Get("beamAnglesDistribution"));
   if (!fAnglesDistribution) {
      LOG(error) << " Beam angular distribution histogram named beamAnglesDistribution not found on the file "
                 << pathToDistributionFile.Data() << "!";
      fBeamAnglesType = DistributionType::kUndefined;
      return;
   }
   fAnglesDistribution->SetDirectory(nullptr);
   LOG(info) << " Beam angular distribution histogram named beamAnglesDistribution found on the file "
             << pathToDistributionFile.Data() << ". Setting of beam angular distribution succesful!";
   fBeamAnglesType = DistributionType::kHistogram;
}

void AtBeamIonGenerator::SetBeamEnergy(double kineticEnergy)
{
   if (kineticEnergy < 0) {
      LOG(error) << " You tried to set the kinetic energy of the beam to a negative value!";
      fKineticEnergyType = DistributionType::kUndefined;
      return;
   }
   fKineticEnergy = kineticEnergy;
   fKineticEnergyType = DistributionType::kDeltaValue;

   LOG(info) << " The kinetic energy of the beam has been set to " << fKineticEnergy << " [MeV/u].";
}

void AtBeamIonGenerator::SetBeamEnergy(double kineticEnergy, double sigmaKE)
{
   SetBeamEnergy(kineticEnergy);

   if (sigmaKE < 0)
      LOG(warning) << " You are setting negative standard deviations for the beam kinetic energy! The simulation code "
                      "should still work fine, but maybe you did not intend to set that value...";
   fSigmaKE = sigmaKE;
   fKineticEnergyType = DistributionType::kGaussian;

   LOG(info) << " Added a standard deviation of the beam kinetic energy of " << fSigmaKE << " [MeV/u].";
}

void AtBeamIonGenerator::SetBeamEnergy(TString pathToDistributionFile)
{
   if (!fFileDistribution)
      fFileDistribution = std::make_unique<TFile>(pathToDistributionFile.Data(), "READ");

   if (fFileDistribution->IsZombie()) {
      LOG(error) << " Tried to pass beam kinetic energy distribution histogram, but the file "
                 << pathToDistributionFile.Data() << " could not be opened!";
      fKineticEnergyType = DistributionType::kUndefined;
      return;
   }

   fKEDistribution = dynamic_cast<TH1F *>(fFileDistribution->Get("kineticEnergyDistribution"));
   if (!fKEDistribution) {
      LOG(error) << " Kinetic energy distribution histogram named kineticEnergyDistribution not found on the file "
                 << pathToDistributionFile.Data() << "!";
      fKineticEnergyType = DistributionType::kUndefined;
      return;
   }
   fKEDistribution->SetDirectory(nullptr);
   LOG(info) << " Kinetic energy distribution histogram named kineticEnergyDistribution found on the file "
             << pathToDistributionFile.Data() << ". Setting of beam angular distribution succesful!";
   fKineticEnergyType = DistributionType::kHistogram;
}

void AtBeamIonGenerator::GenerateBeamIonSourcePoint()
{
   double phi{}, radialDeviation{};
   switch (fBeamProfileType) {
   case DistributionType::kUndefined: LOG(fatal) << " UNDEFINED BEAM PROFILE TYPE! PLEASE FIX ME!"; break;
   case DistributionType::kDeltaValue:
      fOriginX = fSourceX / 10;
      fOriginY = fSourceY / 10;
      fOriginZ = fSourceZ / 10;
      break;
   case DistributionType::kUniform:
      phi = gRandom->Uniform(0, 360) * TMath::DegToRad();
      radialDeviation = gRandom->Uniform(0, std::abs(fBeamRadius)) / 10;

      fOriginX = fSourceX / 10 + radialDeviation * std::cos(phi);
      fOriginY = fSourceY / 10 + radialDeviation * std::sin(phi);
      fOriginZ = fSourceZ / 10;
      break;
   case DistributionType::kGaussian:
      phi = gRandom->Uniform(0, 360) * TMath::DegToRad();
      radialDeviation = gRandom->Gaus(0, std::abs(fBeamRadius)) / 10;

      fOriginX = fSourceX / 10 + radialDeviation * std::cos(phi);
      fOriginY = fSourceY / 10 + radialDeviation * std::sin(phi);
      fOriginZ = fSourceZ / 10;
      break;
   case DistributionType::kHistogram:
      fProfileDistribution->GetRandom2(fOriginX, fOriginY);
      fOriginX /= 10;
      fOriginY /= 10;
      fOriginZ = fSourceZ / 10;
      break;
   }

   LOG(debug) << " Beam source set to (" << fOriginX << ", " << fOriginY << ", " << fOriginZ << ") [cm]...";
}

void AtBeamIonGenerator::GenerateBeamDirection()
{
   switch (fBeamAnglesType) {
   case DistributionType::kUndefined: LOG(fatal) << " UNDEFINED BEAM ANGULAR DISTRIBUTION TYPE! PLEASE FIX ME!"; break;
   case DistributionType::kDeltaValue:
      fThetaXZCurrent = fThetaXZ;
      fThetaYZCurrent = fThetaYZ;
      break;
   case DistributionType::kGaussian:
      fThetaXZCurrent = gRandom->Gaus(fThetaXZ, std::abs(fSigmaXZ));
      fThetaYZCurrent = gRandom->Gaus(fThetaYZ, std::abs(fSigmaYZ));
      break;
   case DistributionType::kHistogram:
      double XZ{}, YZ{};
      fAnglesDistribution->GetRandom2(XZ, YZ);
      fThetaXZCurrent = XZ;
      fThetaYZCurrent = YZ;
      break;
   }

   LOG(debug) << " thetaXZ = " << fThetaXZCurrent << " deg | thetaYZ = " << fThetaYZCurrent << " deg";
}

void AtBeamIonGenerator::GenerateBeamKineticEnergy()
{
   switch (fKineticEnergyType) {
   case DistributionType::kUndefined:
      LOG(fatal) << " UNDEFINED BEAM KINETIC ENERGY DISTRIBUTION TYPE! PLEASE FIX ME!";
      break;
   case DistributionType::kDeltaValue:
      fKineticEnergyCurrent = fKineticEnergy * (fIon->GetMass() + fIon->GetExcEnergy()) * 1000 / fUToMeV;
      break;
   case DistributionType::kGaussian:
      fKineticEnergyCurrent =
         gRandom->Gaus(fKineticEnergy, std::abs(fSigmaKE)) * (fIon->GetMass() + fIon->GetExcEnergy()) * 1000 / fUToMeV;
      break;
   case DistributionType::kHistogram:
      fKineticEnergyCurrent = fKEDistribution->GetRandom() * (fIon->GetMass() + fIon->GetExcEnergy()) * 1000 / fUToMeV;
      break;
   }

   // Technically, the kinetic energy set using these methods could be negative. If that is the case, reset, unless it
   // has undefined kinetic energy distribution type (in which case hopefully crashing the code is a good incentive to
   // make the user fix it).
   if (fKineticEnergyCurrent < 0 && fKineticEnergyType != DistributionType::kUndefined)
      GenerateBeamKineticEnergy();
   else
      LOG(debug) << " Beam kinetic energy set to " << fKineticEnergyCurrent << " [MeV].";
}

void AtBeamIonGenerator::CalculateMomentumVector()
{
   GenerateBeamDirection();
   GenerateBeamKineticEnergy();

   double beamMomentumModulus = std::sqrt(
      fKineticEnergyCurrent / 1000 * (fKineticEnergyCurrent / 1000 + 2 * (fIon->GetMass() + fIon->GetExcEnergy())));

   fPz = beamMomentumModulus / std::sqrt(std::pow(std::tan(fThetaXZCurrent * TMath::DegToRad()), 2) +
                                         std::pow(std::tan(fThetaYZCurrent * TMath::DegToRad()), 2) + 1);

   fPx = fPz * std::tan(fThetaXZCurrent * TMath::DegToRad());
   fPy = fPz * std::tan(fThetaYZCurrent * TMath::DegToRad());

   LOG(debug) << " Beam tri-momentum set to (" << fPx << ", " << fPy << ", " << fPz << ") [GeV/c]";
}

bool AtBeamIonGenerator::ReadEvent(FairPrimaryGenerator *primGen)
{
   bool isBeamEvent = AtVertexPropagator::Instance()->IsBeamEvent();

   if (!isBeamEvent && !fAlwaysSimulateBeam)
      return true;

   if (isBeamEvent) {
      LOG(info) << " Simulating " << fIon->GetName() << " beam ion with PDG code " << fPDGCode << "...";

      GenerateBeamIonSourcePoint();
      CalculateMomentumVector();

      if (fWillSimulateReaction) {
         double beamELoss{};
         if (fMaxELoss < 0)
            beamELoss = gRandom->Uniform(0., fKineticEnergyCurrent);
         else
            beamELoss = gRandom->Uniform(0., fMaxELoss);
         AtVertexPropagator::Instance()->SetRndELoss(beamELoss);
         LOG(info) << cGREEN << " Beam ion energy loss = " << beamELoss << " MeV";
      } else
         AtVertexPropagator::Instance()->SetRndELoss(std::numeric_limits<double>::max());
   }

   primGen->AddTrack(fPDGCode, fPx, fPy, fPz, fOriginX, fOriginY, fOriginZ);

   return true;
}

ClassImp(AtBeamIonGenerator)
