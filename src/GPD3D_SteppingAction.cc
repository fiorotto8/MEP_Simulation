#include "GPD3D_SteppingAction.hh"
#include "GPD3D_EventInfo.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"

GPD3D_SteppingAction::GPD3D_SteppingAction(GPD3D_EventAction* ea, const GPD3D_DetectorConstruction* det)
: G4UserSteppingAction(), fEventAction(ea), fDet(det)
{
}

GPD3D_SteppingAction::~GPD3D_SteppingAction()
{
}

void GPD3D_SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!step) return;

  const G4Track* track = step->GetTrack();
  if (!track) return;

  // We only track containment for the PRIMARY (trackID==1)
  if (track->GetTrackID() != 1) return;

  const G4Event* ev = G4RunManager::GetRunManager()->GetCurrentEvent();
  if (!ev) return;

  auto* info = dynamic_cast<GPD3D_EventInfo*>(ev->GetUserInformation());
  if (!info) return;

  const G4StepPoint* pre  = step->GetPreStepPoint();
  const G4StepPoint* post = step->GetPostStepPoint();
  if (!pre || !post) return;

  const G4VPhysicalVolume* prePV  = pre->GetPhysicalVolume();
  const G4VPhysicalVolume* postPV = post->GetPhysicalVolume();

  const G4LogicalVolume* preLV  = (prePV  ? prePV->GetLogicalVolume()  : nullptr);
  const G4LogicalVolume* postLV = (postPV ? postPV->GetLogicalVolume() : nullptr);

  const bool preInGas  = (preLV  && preLV->GetName()  == "GAS_CHAMBER_LV");
  const bool postInGas = (postLV && postLV->GetName() == "GAS_CHAMBER_LV");

  // Detect entering/exiting gas
  if (!preInGas && postInGas) {
    info->MarkPrimaryEnteredGas();
  }
  if (preInGas && !postInGas) {
    info->MarkPrimaryExitedGas();
  }

  // Detect primary stopping/killed while in gas
  if (preInGas) {
    const auto status = track->GetTrackStatus();
    if (status == fStopAndKill || status == fStopButAlive) {
      info->MarkPrimaryStoppedInGas();
    }
  }

  auto edep = step->GetTotalEnergyDeposit();
  if (edep <= 0) return;

  auto lv = pre->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (!lv) return;

  if (lv == fDet->GetGAGGWall()) {
    fEventAction->edepInGW = true;
    fEventAction->edepGW += edep;
  } else if (lv == fDet->GetGAGGCap()) {
    fEventAction->edepInGC = true;
    fEventAction->edepGC += edep;
  }

}
