/*!
  @file
  @brief Gas cell sensitive detector definition.
*/

#ifndef GPD3D_GASCELLSD_HH
#define GPD3D_GASCELLSD_HH

#include <vector>

#include "G4VSensitiveDetector.hh"

#include "GPD3D_GasCellHit.hh"


class G4Step;
class G4HCofThisEvent;

//! Class representing the gas cell sensitive detector.
/*!
  The hits are accounted in hits in ProcessHits() function which is called
  by Geant4 kernel at each step. A hit is created with each step with non zero
  energy deposit.
 */
class GPD3D_GasCellSD : public G4VSensitiveDetector
{
  public:
    
    //! Constructor.
    GPD3D_GasCellSD(const G4String& name,
                  const G4String& hitsCollectionName);
    
    //! Destructor.
    virtual ~GPD3D_GasCellSD();
    
    //! Overloaded method that creates the hits collection.
    virtual void Initialize(G4HCofThisEvent* hitCollection);
    
    //! Overloaded process-hits method.
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    
    //! Overloaded end-of-event method.
    virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    
    //! The hits collection.
    GPD3D_GasCellHitsCollection* m_hitsCollection;
};

#endif //GPD3D_GASCELLSD_HH
