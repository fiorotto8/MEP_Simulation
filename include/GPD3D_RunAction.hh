#ifndef GPD3D_RUNACTION_HH
#define GPD3D_RUNACTION_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "G4CsvAnalysisManager.hh"
// #include "g4root.hh"

class GPD3D_RunAction : public G4UserRunAction
{
  public:
    GPD3D_RunAction();
    GPD3D_RunAction(const char *);
    virtual ~GPD3D_RunAction();

    // method from the base class
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void SetAnalysis();

  private:
    G4String fName;
};

#endif
