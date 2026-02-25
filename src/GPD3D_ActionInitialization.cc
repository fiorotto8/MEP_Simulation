#include "GPD3D_ActionInitialization.hh"
#include "GPD3D_PrimaryGeneratorAction.hh"
#include "GPD3D_RunAction.hh"
#include "GPD3D_EventAction.hh"


/*!
 */
GPD3D_ActionInitialization::GPD3D_ActionInitialization()
  : G4VUserActionInitialization()
{}


/*!
 */
void GPD3D_ActionInitialization::BuildForMaster() const
{
  SetUserAction(new GPD3D_RunAction());
}


/*!
 */
void GPD3D_ActionInitialization::Build() const
{
  SetUserAction(new GPD3D_PrimaryGeneratorAction());
  SetUserAction(new GPD3D_RunAction());
  SetUserAction(new GPD3D_EventAction());
}
