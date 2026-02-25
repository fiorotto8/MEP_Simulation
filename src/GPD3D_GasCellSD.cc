#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"

#include "GPD3D_GasCellSD.hh"

GPD3D_GasCellSD::GPD3D_GasCellSD(const G4String& name,
                             const G4String& hitsCollectionName) :
  G4VSensitiveDetector(name),
  m_hitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

GPD3D_GasCellSD::~GPD3D_GasCellSD()
{}

void GPD3D_GasCellSD::Initialize(G4HCofThisEvent* hce)
{
  m_hitsCollection = new GPD3D_GasCellHitsCollection(SensitiveDetectorName,
                                                   collectionName[0]);
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, m_hitsCollection);
}

G4bool GPD3D_GasCellSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Total deposit (includes non-ionizing component)
  const G4double eDep = aStep->GetTotalEnergyDeposit();

  // Keep your existing behavior: if gamma and no deposit, skip
  const G4String particle = aStep->GetTrack()->GetDefinition()->GetParticleName();
  if ((eDep == 0.) && (particle == "gamma")) return false;

  // Split into ionizing / non-ionizing (Geant4 provides the non-ion part)
  const G4double eDepNonIon = aStep->GetNonIonizingEnergyDeposit();
  G4double eDepIon = eDep - eDepNonIon;
  if (eDepIon < 0.) eDepIon = 0.; // guard numerical artifacts

  // Kinetic energy at pre-step point
  const G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();

  // Metadata
  const G4int pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  const G4double time = aStep->GetPreStepPoint()->GetGlobalTime();

  // Process type/subtype (numeric, avoid strings in ROOT)
  G4int creatorType = -1;
  G4int creatorSubType = -1;
  if (aStep->GetTrack()->GetCreatorProcess() != nullptr) {
    creatorType = aStep->GetTrack()->GetCreatorProcess()->GetProcessType();
    creatorSubType = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();
  }

  G4int stepType = -1;
  G4int stepSubType = -1;
  const G4VProcess* stepProc = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if (stepProc != nullptr) {
    stepType = stepProc->GetProcessType();
    stepSubType = stepProc->GetProcessSubType();
  }

  GPD3D_GasCellHit* newHit = new GPD3D_GasCellHit();

  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetParentID(aStep->GetTrack()->GetParentID());
  newHit->SetStepNum(aStep->GetTrack()->GetCurrentStepNumber());

  newHit->SetEdep(eDep);
  newHit->SetEdepIon(eDepIon);
  newHit->SetEdepNonIon(eDepNonIon);

  newHit->SetEnergy(energy);

  newHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
  newHit->SetMomDir(aStep->GetPreStepPoint()->GetMomentumDirection());

  // Keep existing strings (used elsewhere). ROOT output will not use them.
  newHit->SetTrackProcName(aStep->GetTrack()->GetCreatorModelName());
  if (stepProc != nullptr) {
    newHit->SetStepProcName(stepProc->GetProcessName());
  } else {
    newHit->SetStepProcName("");
  }
  newHit->SetPartName(particle);

  newHit->SetPDG(pdg);
  newHit->SetGlobalTime(time);

  newHit->SetCreatorProcType(creatorType);
  newHit->SetCreatorProcSubType(creatorSubType);
  newHit->SetStepProcType(stepType);
  newHit->SetStepProcSubType(stepSubType);

  newHit->SetTrkLength(aStep->GetTrack()->GetTrackLength());
  newHit->SetPartState(aStep->GetTrack()->GetTrackStatus());

  // Ionization track sampling still left as-is (commented in your code)
  // newHit->SetIonTrack(GPD3D_DetectorSimulationSvc::sampleInStep(
  //                     aStep, 1./GPD3D_RunConfiguration::fanoFactor,
  //                     GPD3D_RunConfiguration::ionizationEnergy));

  newHit->SetFirstStep(aStep->IsFirstStepInVolume());
  newHit->SetLastStep(aStep->IsLastStepInVolume());

  newHit->SetVertPos(aStep->GetTrack()->GetVertexPosition());
  newHit->SetVertDir(aStep->GetTrack()->GetVertexMomentumDirection());
  newHit->SetVertEnergy(aStep->GetTrack()->GetVertexKineticEnergy());

  m_hitsCollection->insert(newHit);
  return true;
}

void GPD3D_GasCellSD::EndOfEvent(G4HCofThisEvent*)
{
}
