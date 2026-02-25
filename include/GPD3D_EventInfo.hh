#ifndef GPD3D_EVENTINFO_HH
#define GPD3D_EVENTINFO_HH

#include "G4VUserEventInformation.hh"
#include "G4ios.hh"

/**
 * Lightweight per-event flags (optional).
 * You can extend/actually set these flags in SteppingAction/SD if you want.
 */
class GPD3D_EventInfo : public G4VUserEventInformation {
public:
  GPD3D_EventInfo() { Reset(); }
  ~GPD3D_EventInfo() override = default;

  // Convenience setters (used by SteppingAction)
  void MarkPrimaryEnteredGas()   { primaryEnteredGas   = true; }
  void MarkPrimaryExitedGas()    { primaryExitedGas    = true; }
  void MarkPrimaryStoppedInGas() { primaryStoppedInGas = true; }

  void Reset() {
    primaryEnteredGas   = false;
    primaryExitedGas    = false;
    primaryStoppedInGas = false;
  }

  bool PrimaryContainedInGas() const {
    // “contained” meaning: entered + stopped and did not exit
    return primaryEnteredGas && primaryStoppedInGas && !primaryExitedGas;
  }

  void Print() const override {
    G4cout << "GPD3D_EventInfo: "
           << "enteredGas=" << primaryEnteredGas
           << " exitedGas=" << primaryExitedGas
           << " stoppedInGas=" << primaryStoppedInGas
           << " containedGas=" << PrimaryContainedInGas()
           << G4endl;
  }

  bool primaryEnteredGas;
  bool primaryExitedGas;
  bool primaryStoppedInGas;
};

#endif
