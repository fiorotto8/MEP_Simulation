#ifndef GPD3D_PHYSICS_LIST_HH
#define GPD3D_PHYSICS_LIST_HH

#include <G4VModularPhysicsList.hh>
#include "globals.hh"

class G4VPhysicsConstructor;

class GPD3D_PhysicsList : public G4VModularPhysicsList
{
public:
  GPD3D_PhysicsList();
  ~GPD3D_PhysicsList();

  //! Optional virtual methods, to gain direct control on 
  //! the particle/processes definition. Not used here
  /*
  void 	ConstructParticle () override;
  void 	ConstructProcess () override;
  */

  //! Mandatory method 
  void 	SetCuts ();

};

#endif
