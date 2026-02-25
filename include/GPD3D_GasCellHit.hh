/*!
  @file
  @brief Gas cell hit definition.
*/

#ifndef GPD3D_GASCELLHIT_HH
#define GPD3D_GASCELLHIT_HH

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
#include "tls.hh"

#include "GPD3D_IonizationTrack.hh"

//! Gas cell hit class.
class GPD3D_GasCellHit : public G4VHit
{
  public:

    GPD3D_GasCellHit();
    GPD3D_GasCellHit(const GPD3D_GasCellHit&);
    virtual ~GPD3D_GasCellHit();

    const GPD3D_GasCellHit& operator=(const GPD3D_GasCellHit&);
    G4int operator==(const GPD3D_GasCellHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // ---- existing setters ----
    void SetTrackID (G4int track) {m_trackID = track;}
    void SetParentID (G4int parent) {m_parentID = parent;}
    void SetStepNum (G4int step) {m_stepNum = step;}
    void SetEdep (G4double de) {m_eDep = de;}
    void SetEnergy (G4double energy) {m_energy = energy;}
    void SetTrkLength (G4double length) {m_trkLength = length;}
    void SetPos (G4ThreeVector xyz) {m_pos = xyz;}
    void SetMomDir (G4ThreeVector xyz) {m_momDir = xyz;}
    void SetTrackProcName (G4String proc) {m_trackProcName = proc;}
    void SetStepProcName (G4String proc) {m_stepProcName = proc;}
    void SetPartName (G4String part) {m_partName = part;}
    void SetPartState (G4TrackStatus state) {m_partState = state;}
    void SetIonTrack (GPD3D_IonizationTrack ion) {m_ionTrack = ion;}
    void SetFirstStep (G4bool flag) {m_firstStep = flag;}
    void SetLastStep (G4bool flag) {m_lastStep = flag;}
    void SetVertPos (G4ThreeVector xyz) {m_vertPos = xyz;}
    void SetVertDir (G4ThreeVector xyz) {m_vertDir = xyz;}
    void SetVertEnergy (G4double energy) {m_vertEnergy = energy;}

    // ---- NEW numeric-only metadata (for ROOT, avoid strings) ----
    void SetEdepIon(G4double de) { m_eDepIon = de; }
    void SetEdepNonIon(G4double de) { m_eDepNonIon = de; }
    void SetPDG(G4int pdg) { m_pdg = pdg; }
    void SetGlobalTime(G4double t) { m_globalTime = t; }

    void SetCreatorProcType(G4int t) { m_creatorProcType = t; }
    void SetCreatorProcSubType(G4int t) { m_creatorProcSubType = t; }
    void SetStepProcType(G4int t) { m_stepProcType = t; }
    void SetStepProcSubType(G4int t) { m_stepProcSubType = t; }

    // ---- existing getters ----
    G4int GetTrackID() const {return m_trackID;}
    G4int GetParentID() const {return m_parentID;}
    G4int GetStepNum() const {return m_stepNum;}
    G4double GetEdep() const {return m_eDep;}
    G4double GetEnergy() const {return m_energy;}
    G4double GetTrkLength() const {return m_trkLength;}
    G4ThreeVector GetPos() const {return m_pos;}
    G4ThreeVector GetMomDir() const {return m_momDir;}
    G4String GetTrackProcName() const {return m_trackProcName;}
    G4String GetStepProcName() const {return m_stepProcName;}
    G4String GetPartName() const {return m_partName;}
    G4TrackStatus GetPartState() const {return m_partState;}
    GPD3D_IonizationTrack GetIonTrack() const {return m_ionTrack;}
    G4bool IsFirstStep() const {return m_firstStep;}
    G4bool IsLastStep() const {return m_lastStep;}
    G4ThreeVector GetVertPos() const {return m_vertPos;}
    G4ThreeVector GetVertDir() const {return m_vertDir;}
    G4double GetVertEnergy() const {return m_vertEnergy;}

    // ---- NEW getters ----
    G4double GetEdepIon() const { return m_eDepIon; }
    G4double GetEdepNonIon() const { return m_eDepNonIon; }
    G4int GetPDG() const { return m_pdg; }
    G4double GetGlobalTime() const { return m_globalTime; }

    G4int GetCreatorProcType() const { return m_creatorProcType; }
    G4int GetCreatorProcSubType() const { return m_creatorProcSubType; }
    G4int GetStepProcType() const { return m_stepProcType; }
    G4int GetStepProcSubType() const { return m_stepProcSubType; }

  private:

    G4int m_trackID;
    G4int m_parentID;
    G4int m_stepNum;
    G4double m_eDep;
    G4double m_energy;
    G4double m_trkLength;
    G4ThreeVector m_pos;
    G4ThreeVector m_momDir;
    G4String m_trackProcName;
    G4String m_stepProcName;
    G4String m_partName;
    G4TrackStatus m_partState;
    GPD3D_IonizationTrack m_ionTrack;
    G4bool m_firstStep;
    G4bool m_lastStep;

    G4ThreeVector m_vertPos;
    G4ThreeVector m_vertDir;
    G4double m_vertEnergy;

    // ---- NEW numeric-only metadata ----
    G4double m_eDepIon;
    G4double m_eDepNonIon;
    G4int    m_pdg;
    G4double m_globalTime;

    G4int    m_creatorProcType;
    G4int    m_creatorProcSubType;
    G4int    m_stepProcType;
    G4int    m_stepProcSubType;
};

typedef G4THitsCollection<GPD3D_GasCellHit> GPD3D_GasCellHitsCollection;

extern G4ThreadLocal G4Allocator<GPD3D_GasCellHit>* GPD3D_GasCellHitAllocator;

inline void* GPD3D_GasCellHit::operator new(size_t)
{
  if(!GPD3D_GasCellHitAllocator)
      GPD3D_GasCellHitAllocator = new G4Allocator<GPD3D_GasCellHit>;
  return (void *) GPD3D_GasCellHitAllocator->MallocSingle();
}

inline void GPD3D_GasCellHit::operator delete(void *hit)
{
  GPD3D_GasCellHitAllocator->FreeSingle((GPD3D_GasCellHit*) hit);
}

#endif //GPD3D_GASCELLHIT_HH
