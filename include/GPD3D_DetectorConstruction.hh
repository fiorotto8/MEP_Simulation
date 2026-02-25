#ifndef GPD3D_DETECTORCONSTRUCTION_HH
#define GPD3D_DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UniformElectricField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "globals.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class G4VPhysicalVolume;
class G4LogicalVolume;

class GPD3D_DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    GPD3D_DetectorConstruction();

    virtual ~GPD3D_DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    //! Overloaded method for sensitive detectors and fields construction.
    virtual void ConstructSDandField();

    // Methods to retriev pointer to GAGG wall/cap logical volumes, for use in event action (to distinguish wall vs cap hits)
    G4LogicalVolume* GetGAGGWall() const { return lvGAGGWall; }
    G4LogicalVolume* GetGAGGCap() const { return lvGAGGCap; }
    
  /*
  private:
    //! Construct the electric field.
    void ElectricFieldSetup();

    //! The gas cell logical volume.
    G4LogicalVolume* gas_chamber_Log;
    
    //! Electric field.
    G4ElectricField* m_EMfield;
    
    //! The electric field equation of motion.
    G4EqMagElectricField* m_equation;
    
    //! Stepper for equation of motion integration.
    G4MagIntegratorStepper* m_stepper;
    
    //! Minimum step size used during integration.
    G4double m_minStep;
    
    //! Chord finder.
    G4ChordFinder* m_chordFinder;
  */
  private:
    void ElectricFieldSetup();

    G4LogicalVolume* gas_chamber_Log = nullptr;
    G4LogicalVolume* gas_sensitive_Log = nullptr;

    G4LogicalVolume* lvGAGGWall = nullptr;
    G4LogicalVolume* lvGAGGCap = nullptr;


    G4ElectricField* m_EMfield = nullptr;
    G4EqMagElectricField* m_equation = nullptr;
    G4MagIntegratorStepper* m_stepper = nullptr;
    G4double m_minStep = 0.0;
    G4ChordFinder* m_chordFinder = nullptr;

    
};


#endif

 