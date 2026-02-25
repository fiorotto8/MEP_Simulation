/*!
  @file
  @brief Physical process and particles definition.
*/

#ifndef GPD3D_PHYSICS_H
#define GPD3D_PHYSICS_H

#include "G4VPhysicsConstructor.hh"


//! Class describing the physical process and particles for simulation.
/*!
  This is essentially taken from the G4EmLivermorePolarizedPhysics constructor:
  http://geant4.cern.ch/support/source/geant4/source/physics_lists/constructors/electromagnetic/src/G4EmLivermorePolarizedPhysics.cc
 */
class GPD3D_Physics : public G4VPhysicsConstructor
{
  public:
    
    //! Constructor.
    explicit GPD3D_Physics(G4int ver=1, const G4String& name="");
    
    //! Destructor.
    virtual ~GPD3D_Physics();
    
    //! Overloaded costruct-particle method.
    virtual void ConstructParticle();
    
    //! Overloaded costruct-process method.
    virtual void ConstructProcess();

  private:
    
    //! Verbose level
    G4int  verbose;
  
};

#endif //GPD3D_PHYSICS_H
