#ifndef GPD3D_STEPPINGACTION_HH
#define GPD3D_STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "GPD3D_DetectorConstruction.hh"
#include "GPD3D_EventAction.hh"

class GPD3D_SteppingAction : public G4UserSteppingAction
{
public:
  GPD3D_SteppingAction(GPD3D_EventAction* ea, const GPD3D_DetectorConstruction* det);
  virtual ~GPD3D_SteppingAction();

  virtual void UserSteppingAction(const G4Step* step) override;

private:
  GPD3D_EventAction* fEventAction;
  const GPD3D_DetectorConstruction* fDet;

};

#endif
