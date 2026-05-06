/*
   I have been having issues with the generator classes that already exist, and in my opinion their logic is too
   obscure. Therefore, I have decided to create my own simulation classes so I fully understand their inner workings...
*/

#ifndef ATBEAMIONGENERATOR_H
#define ATBEAMIONGENERATOR_H

#include <FairGenerator.h>

#include <Rtypes.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>

#include <memory>

class FairPrimaryGenerator;
class FairIon;
class TBuffer;
class TClass;
class TMemberInspector;

class AtBeamIonGenerator : public FairGenerator {
public:
   // Enumerating different cases for the (X,Y) profile of the beam, its angular distribution and energy distribution.
   // Here, we are assuming that there is no correlation in (X,Y) position of ion, direction and energy, which is most
   // likely false... TO-DO: Implement the case where there is correlation in position, beam direction and energy.
   enum class DistributionType { kUndefined = 0, kDeltaValue = 1, kGaussian = 2, kUniform = 3, kHistogram = 4 };

protected:
   // TFile containing distribution histograms.
   std::unique_ptr<TFile> fFileDistribution{nullptr};

   // Beam profile information.
   DistributionType fBeamProfileType{DistributionType::kDeltaValue};
   double fSourceX{}, fSourceY{}, fSourceZ{}; // [cm] x source, y source, z source.
   double fBeamRadius{};                      // [cm] Beam profile radius, which may be used in case kGaussian.
   TH2F *fProfileDistribution{nullptr}; // 2D histogram of the beam (X,Y) profile, both in [cm].

   // Beam direction information.
   DistributionType fBeamAnglesType{DistributionType::kDeltaValue};
   double fThetaXZ{}, fThetaYZ{};      // [deg] Angles between beam and Z direction.
   double fSigmaXZ{}, fSigmaYZ{};      // [deg] Angles distributions standard deviations, which may be used in
                                       // case kGaussian.
   TH2F *fAnglesDistribution{nullptr}; // 2D histogram of the beam angular distribution
                                                       // (X-Z and Y-Z angles), both in [deg].

   // Beam energy information.
   DistributionType fKineticEnergyType{DistributionType::kUndefined};
   double fKineticEnergy{-1};      // [MeV/u] to be used in case kDeltaValue or kGaussian.
   double fSigmaKE{};              // [MeV/u] Kinetic energy standard deviation, which may be used in case kGaussian.
   TH1F *fKEDistribution{nullptr}; // 1D histogram of the kinetic energy distribution, in [MeV/u].

   // Particle information.
   FairIon *fIon;                        // Pointer to the FairIon to be generated.
   int fZ{-1}, fA{-1}, fChargeState{-1}; // Atomic number, mass number and charge state of the ion to be generated.
   double fMass{-1};                     // [MeV] Ground state mass of ion.
   double fEx{};                         // [MeV] Excitation energy of ion.
   TString fIonName{"undefined"};
   int fPDGCode{-99999};

   // Event by event information.
   double fOriginX{-999}, fOriginY{-999}, fOriginZ{-999}; // [cm] Origin of ion.
   double fThetaXZCurrent{}, fThetaYZCurrent{};           // [deg] Direction angles.
   double fKineticEnergyCurrent{-1};                      // [MeV] KE of ion.
   double fPx{-999}, fPy{-999}, fPz{-999};                // [MeV] Momentum of ion.

   // Adding ion logic and vertex propagator configuration.
   bool IonAlreadyExists(TString ionName);
   bool AddIonToFairRunSim();
   void ConfigureVertexPropagator();

   // Private beam ion generation methods
   void GenerateBeamIonSourcePoint();
   void GenerateBeamDirection();
   void GenerateBeamKineticEnergy();
   void CalculateMomentumVector();

   // u to MeV conversion constant.
   double fUToMeV{931.4936};

   // Other configuration parameters.
   bool fAlwaysSimulateBeam{false};
   bool fWillSimulateReaction{true};
   double fMaxELoss{-1}; // [MeV] Maximum energy loss before reation happens.

public:
   /** Constructor using beam ion information only. Optionally, you may pass kinetic energy value and its standard
    *deviation. For histogram distribution of kinetic energy values, use the SetBeamEnergy() method.
    **@param ionName  Ion name
    **@param pdgCode PDG Code of the ion
    **@param Z  Atomic number
    **@param A  Mass number
    **@param chargeState  Charge state number
    **@param mass  G.S. mass of ion in [MeV]
    **@param Ex  Excitation energy of ion in [MeV] (0 by default)
    **@param kineticEnergy  Kinetic energy of ion in [MeV/u] (0 by default)
    **@param sigmaKE  Standard deviation of excitation energy of ion in [MeV/u] (0 by default)
    **/
   AtBeamIonGenerator(TString ionName, int pdgCode, int Z, int A, int chargeState, double mass, double Ex = 0,
                      double kineticEnergy = -1, double sigmaKE = -1);

   /** Constructor using beam ion name only. Optionally, you may pass kinetic energy value and its standard
    *deviation. For histogram distribution of kinetic energy values, use the SetBeamEnergy() method.
    **@param ionName  Ion name
    **@param pdgCode PDG Code of the ion
    **@param kineticEnergy  Kinetic energy of ion in [MeV/u] (0 by default)
    **@param sigmaKE  Standard deviation of excitation energy of ion in [MeV/u] (0 by default)
    **/
   AtBeamIonGenerator(TString ionName, int pdgCode, double kineticEnergy = -1, double sigmaKE = -1);

   /** Destructor **/
   virtual ~AtBeamIonGenerator();

   /** Modifiers **/
   void SetBeamSource(double z, double x, double y);
   void SetBeamSource(double z, double x, double y, double beamRadius);
   void SetBeamSourceGaussian(double z, double x, double y, double beamRadius);
   void SetBeamSource(double z, TString pathToDistributionFile);

   void SetBeamDirection(double thetaXZ, double thetaYZ);
   void SetBeamDirection(double thetaXZ, double thetaYZ, double sigmaXZ, double sigmaYZ);
   void SetBeamDirection(TString pathToDistributionFile);

   void SetBeamEnergy(double kineticEnergy);
   void SetBeamEnergy(double kineticEnergy, double sigmaKE);
   void SetBeamEnergy(TString pathToDistributionFile);

   void SetAlwaysSimulateBeam(bool value = true) { fAlwaysSimulateBeam = value; }
   void SetWillSimulateReaction(bool value = true) { fWillSimulateReaction = value; }
   void SetMaxEnergyLoss(double value) { fMaxELoss = value; }

   /** Method ReadEvent
   ** Generates the specified ion and hands it to the FairPrimaryGenerator.
   **/
   virtual bool ReadEvent(FairPrimaryGenerator *primGen) override;

   ClassDefOverride(AtBeamIonGenerator, 1)
};

#endif
