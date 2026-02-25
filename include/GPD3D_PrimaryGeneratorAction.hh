#ifndef GPD3D_PRIMARYGENERATORACTION_HH
#define GPD3D_PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "G4GenericMessenger.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"

#include <vector>

class GPD3D_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    // Phase C: sphere injection for 4pi isotropic environment
    G4bool m_useSphere;
    G4double m_sphereRadius;     // length units (e.g. mm)
    G4ThreeVector m_sphereCenter;
    // Beam-mode controls (used when useSphere == false)
    G4double      m_beamSigma;    // Gaussian sigma in length units
    G4ThreeVector m_beamCenter;   // center of Gaussian in XY (Z handled separately)
    G4double      m_beamZ;        // fixed Z for beam source
    G4bool m_useGpsMacro;
    G4double m_posThetaMin; // radians
    G4double m_posThetaMax; // radians
    
    G4GenericMessenger* m_genMessenger;

    GPD3D_PrimaryGeneratorAction();    
    virtual ~GPD3D_PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);      

  private:
    //! Convert theta and phi into a vector.
    static G4ThreeVector directionVector(double theta, double phi);

    //! Convert a polarization angle to a three-vector.
    static G4ThreeVector polarizationVector(double polAngle);
    //! Setup the polarization of the source.
    // static void setupPolarization(G4GeneralParticleSource *source);

    static G4ThreeVector direction;

    // --- Spectrum energy (CSV) ---
    void SetSpectrumFile(const G4String& path);
    void LoadSpectrumFromCsv(const G4String& path);
    G4double SampleSpectrumEnergy() const;
    
    G4bool   m_useSpectrum;
    G4String m_spectrumFile;

    // Spectrum tables:
    // m_specE     : energies (in internal Geant4 units)
    // m_specW     : weights at points
    // m_specSegA  : per-segment exponent a for log-log interpolation
    // m_specCum   : normalized CDF at segment boundaries (size = nseg+1)
    std::vector<G4double> m_specE;
    std::vector<G4double> m_specW;
    std::vector<G4double> m_specSegA;
    std::vector<G4double> m_specCum;
    
    G4GeneralParticleSource* fGPS;
};

#endif
