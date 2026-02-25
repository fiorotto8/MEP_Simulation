// absolute z extension of gas volume 
//z_min = 28.590 - 14.997 = 13.593 mm
//z_max = 28.590 + 14.997 = 43.587 mm
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4ClassicalRK4.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UniformElectricField.hh"
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
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4Transform3D.hh"
#include "G4DisplacedSolid.hh"


#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4SDManager.hh>
// #include <G4Colour.hh>

#include <sstream>
#include <math.h>

using namespace std;

#include "GPD3D_DetectorConstruction.hh"
#include "GPD3D_GasCellSD.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"

GPD3D_DetectorConstruction::GPD3D_DetectorConstruction()
: G4VUserDetectorConstruction()
{
}

GPD3D_DetectorConstruction::~GPD3D_DetectorConstruction()
{
}

G4VPhysicalVolume* GPD3D_DetectorConstruction::Construct()
{
  auto nist = G4NistManager::Instance();

  // ---------- WORLD ----------
  const G4double worldSize = 1.0*m;

  auto worldBox = new G4Box("world", worldSize/2, worldSize/2, worldSize/2);

  auto vacuum = new G4Material("Vacuum", 1., 1.01*g/mole,
                               universe_mean_density, kStateGas,
                               0.1*kelvin, 1.e-19*pascal);

  auto worldLV = new G4LogicalVolume(worldBox, vacuum, "world");
  auto worldVis = new G4VisAttributes();
  worldVis->SetVisibility(false);
  worldLV->SetVisAttributes(worldVis);

  auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0);

  // ---------- MATERIALS ----------
  auto H  = nist->FindOrBuildElement("H");
  auto C  = nist->FindOrBuildElement("C");
  auto O  = nist->FindOrBuildElement("O");
  auto ArE = nist->FindOrBuildElement("Ar");
  auto Si = nist->FindOrBuildMaterial("G4_Si");
  auto Al = nist->FindOrBuildMaterial("G4_Al");
  auto Ti = nist->FindOrBuildMaterial("G4_Ti");
  auto Cu = nist->FindOrBuildMaterial("G4_Cu");
  auto Be = nist->FindOrBuildMaterial("G4_Be");
  auto Kapton = nist->FindOrBuildMaterial("G4_KAPTON");

  const G4double Tgas = 293.15*kelvin;
  const G4double Pgas = 3.0*atmosphere;

  // DME
  const G4double dmeDensity = 3.0 * 1.9e-3 * g/cm3;
  auto dme = new G4Material("DME", dmeDensity, 3, kStateGas, Tgas, Pgas);
  dme->AddElement(C, 2);
  dme->AddElement(H, 6);
  dme->AddElement(O, 1);
  dme->GetIonisation()->SetMeanExcitationEnergy(60.*eV);

  // Argon
  const G4double ArDensity = 3.0 * 1.7836 * mg/cm3;
  auto argon = new G4Material("Argon", ArDensity, 1, kStateGas, Tgas, Pgas);
  argon->AddElement(ArE, 1);

  // Mixture
  const G4double mixDensity = 3.0 * 1.809 * mg/cm3;
  auto mix = new G4Material("Mixture", mixDensity, 2, kStateGas, Tgas, Pgas);
  mix->AddMaterial(dme, 0.22);
  mix->AddMaterial(argon, 0.78);

  // PEEK
  auto PEEK = new G4Material("PEEK", 1.32*g/cm3, 3);
  PEEK->AddElement(C, 19);
  PEEK->AddElement(H, 12);
  PEEK->AddElement(O, 3);

  // FR4 (simple)
  auto Epoxy = new G4Material("Epoxy", 1.20*g/cm3, 3);
  Epoxy->AddElement(C, 0.70);
  Epoxy->AddElement(H, 0.06);
  Epoxy->AddElement(O, 0.24);

  auto Glass = new G4Material("GlassFibre", 2.60*g/cm3, 2);
  Glass->AddElement(nist->FindOrBuildElement("Si"), 1);
  Glass->AddElement(O, 2);

  auto FR4 = new G4Material("FR4", 1.85*g/cm3, 2);
  FR4->AddMaterial(Epoxy, 0.38);
  FR4->AddMaterial(Glass, 0.62);

  // SU8 pillars
  auto SU8 = new G4Material("SU8", 1.2*g/cm3, 3);
  SU8->AddElement(C, 87);
  SU8->AddElement(H, 118);
  SU8->AddElement(O, 16);

  //GAGG 
  auto Gd  = nist->FindOrBuildElement("Gd");
  auto AlE = nist->FindOrBuildElement("Al");
  auto Ga  = nist->FindOrBuildElement("Ga");

  const G4double GAGG_density = 6.63*g/cm3; //:contentReference[oaicite:2]{index=2}
  auto GAGG = new G4Material("GAGG", GAGG_density, 4);
  GAGG->AddElement(Gd, 3);
  GAGG->AddElement(AlE, 2);
  GAGG->AddElement(Ga, 3);
  GAGG->AddElement(O, 12);


  // ---------- VIS ----------
  auto visGrey  = new G4VisAttributes(G4Colour::Grey());  visGrey->SetForceSolid(true);
  auto visRed   = new G4VisAttributes(G4Colour::Red());   visRed->SetForceSolid(true);
  auto visBlue  = new G4VisAttributes(G4Colour::Blue());  visBlue->SetForceSolid(true);
  auto visGreen = new G4VisAttributes(G4Colour(0.4,0.9,0.0,0.4)); visGreen->SetForceSolid(true);
  auto visGold  = new G4VisAttributes(G4Colour(1.0,0.9,0.1,0.4)); visGold->SetForceSolid(true);
  auto visWhite = new G4VisAttributes(G4Colour::White()); visWhite->SetForceSolid(true);
  auto visBrown = new G4VisAttributes(G4Colour(0.6,0.3,0.1,0.6)); visBrown->SetForceSolid(true);
  auto visSens = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.25)); // transparent red
  visSens->SetForceSolid(true);


  // ---------- GEOMETRY (all along +Z, no overlaps by construction) ----------
  const G4double gap = 5.0*um;          // separation between touching faces
  const G4double sub_margin = 1.0*um;   // boolean margins

  // Titanium frame (disk - square hole)
  const G4double Ti_aperture_square_side = 16.0*mm;
  const G4double Ti_frame_disk_dia   = 35.0*2.0*mm;
  const G4double Ti_frame_disk_thick = 5.0*mm;

  // Be window (square)
  const G4double Be_window_square_side = 22.0*mm;
  const G4double Be_window_thick       = 200.0*um;

  // PEEK base + groove
  const G4double PEEK_outer_dia        = 54.0*2.0*mm;
  const G4double PEEK_groove_inner_dia = Ti_frame_disk_dia;
  const G4double PEEK_groove_thick     = Ti_frame_disk_thick;
  const G4double PEEK_base_inner_dia   = 20.0*2.0*mm;
  const G4double PEEK_base_thick       = 30.0*mm;

  // Cu anode (disk - square hole)
  const G4double Cu_Anode_disk_dia     = 30.0*2.0*mm;
  const G4double Cu_Anode_disk_thick   = 500.0*um;
  const G4double Cu_Anode_aperture_side= 15.0*mm;

  // “Vespel” (Kapton) base + groove
  const G4double Vespel_outer_dia        = PEEK_outer_dia;
  const G4double Vespel_groove_inner_dia = Cu_Anode_disk_dia;
  const G4double Vespel_groove_thick     = Cu_Anode_disk_thick;
  const G4double Vespel_base_thick       = 5.5*mm;
  const G4double Vespel_base_ap_side     = 20.0*mm;

  // PCB + Al support
  const G4double PCB_dia = Vespel_outer_dia;
  const G4double PCB_thick = 1.57*mm;
  const G4double PCBSupportAlThick = 6.0*mm;

  // Timepix / foil / pillars (your values)
  const G4double thickness_timepix = 0.013*mm;
  const G4double foilThickness     = 0.001*mm;
  const G4double foilSizeXY        = 14.08*mm;
  const G4double pillarDiameter    = 30.0*um;
  const G4double pillarHeight      = 50.0*um;
  const G4double pillarSpacing     = 110.0*um;
  const G4int numPillarsX = static_cast<G4int>(foilSizeXY / pillarSpacing);
  const G4int numPillarsY = static_cast<G4int>(foilSizeXY / pillarSpacing);
  
  // ---------- GAGG SHIELD PARAMETERS ----------
  const G4double shieldThick = 10.0*mm;       // 1 cm thick
  const G4double shieldClear = 0.5*mm;        // small clearance from detector
  const G4double detRmax     = PEEK_outer_dia/2.;  // largest radius in your detector
  const G4double shieldRin   = detRmax + shieldClear;
  const G4double shieldRout  = shieldRin + shieldThick;

  // ---------- AL BOX SHIELD PARAMETERS ----------
  const G4double alBoxOuterXY = 15.0*cm;     // outer size in X and Y
  const G4double alBoxOuterZ  = 15.0*cm;     // outer size in Z
  const G4double alBoxThick   = 1.5*mm;      // wall thickness

  // hole on +Z face (choose one)
  const G4double alHoleRadius = 1.0*cm;      // example: 3 cm radius
  // if you prefer a square hole:
  // const G4double alHoleSide   = 6.0*cm;


  // ---------- SOLIDS ----------
  auto sPCBSupport = new G4Tubs("PCBSupportSolid", 0., PCB_dia/2., PCBSupportAlThick/2., 0.*deg, 360.*deg);
  auto sPCB        = new G4Tubs("PCBSolid",        0., PCB_dia/2., PCB_thick/2.,        0.*deg, 360.*deg);

  // Vespel base: disk - square hole
  auto sVespelBaseOuter = new G4Tubs("OuterVespelBase", 0., Vespel_outer_dia/2., Vespel_base_thick/2., 0.*deg, 360.*deg);
  auto sVespelBaseInner = new G4Box("InnerVespelBase", Vespel_base_ap_side/2., Vespel_base_ap_side/2., Vespel_base_thick/2. + sub_margin);
  auto sVespelBase = new G4SubtractionSolid("VespelBaseSolid", sVespelBaseOuter, sVespelBaseInner);

  // Vespel groove ring (radially disjoint from anode disk)
  auto sVespelGroove = new G4Tubs("VespelGrooveSolid", Vespel_groove_inner_dia/2., Vespel_outer_dia/2., Vespel_groove_thick/2., 0.*deg, 360.*deg);

  // Cu anode: disk - square hole
  auto sCuAnodeOuter = new G4Tubs("OuterCuAnode", 0., Cu_Anode_disk_dia/2., Cu_Anode_disk_thick/2., 0.*deg, 360.*deg);
  auto sCuAnodeInner = new G4Box("InnerCuAnode", Cu_Anode_aperture_side/2., Cu_Anode_aperture_side/2., Cu_Anode_disk_thick/2. + sub_margin);
  auto sCuAnode = new G4SubtractionSolid("CuAnode", sCuAnodeOuter, sCuAnodeInner);

  // PEEK base ring
  auto sPEEKBase = new G4Tubs("PEEKBaseSolid", PEEK_base_inner_dia/2., PEEK_outer_dia/2., PEEK_base_thick/2., 0.*deg, 360.*deg);

  // PEEK groove ring (radially disjoint from Ti frame disk)
  auto sPEEKGroove = new G4Tubs("PEEKGrooveSolid", PEEK_groove_inner_dia/2., PEEK_outer_dia/2., PEEK_groove_thick/2., 0.*deg, 360.*deg);

  // Ti frame: disk - square hole
  auto sTiOuter = new G4Tubs("OuterTiFrame", 0., Ti_frame_disk_dia/2., Ti_frame_disk_thick/2., 0.*deg, 360.*deg);
  auto sTiInner = new G4Box("InnerBox", Ti_aperture_square_side/2., Ti_aperture_square_side/2., Ti_frame_disk_thick/2. + sub_margin);
  auto sTiFrame = new G4SubtractionSolid("TiFrame", sTiOuter, sTiInner);

  // Be window (square plate)
  auto sBeWin = new G4Box("BeWindowSolid", Be_window_square_side/2., Be_window_square_side/2., Be_window_thick/2.);

  // Timepix / foil
  auto sTimepix = new G4Box("grid_pix", foilSizeXY/2., foilSizeXY/2., thickness_timepix/2.);
  auto sFoil    = new G4Box("foil_al",  foilSizeXY/2., foilSizeXY/2., foilThickness/2.);

  // GAS volume: inner cylinder of PEEK base (slightly smaller to avoid boundary warnings)
  const G4double gasR = (PEEK_base_inner_dia/2.) - 2.0*um;
  auto sGas = new G4Tubs("GasDriftSolid", 0., gasR, PEEK_base_thick/2. - 2.0*um, 0.*deg, 360.*deg);

  // Sensitive gas volume: a box above the Timepix area, spanning full drift thickness
  const G4double sensMarginXY = 1.0*um;  // small margin to avoid boundary effects
  const G4double sensMarginZ  = 1.0*um;

  auto sGasSens = new G4Box("GasSensitiveSolid",
                            foilSizeXY/2. - sensMarginXY,
                            foilSizeXY/2. - sensMarginXY,
                            (PEEK_base_thick/2. - 2.0*um) - sensMarginZ);


  // ---------- LOGICAL VOLUMES ----------
  auto lvPCBSupport = new G4LogicalVolume(sPCBSupport, Al, "PCBSupportAl_LV");
  auto lvPCB        = new G4LogicalVolume(sPCB, FR4, "PCB_LV");
  auto lvVespelBase = new G4LogicalVolume(sVespelBase, Kapton, "VespelBase_LV");
  auto lvVespelGro  = new G4LogicalVolume(sVespelGroove, Kapton, "VespelGroove_LV");
  auto lvCuAnode    = new G4LogicalVolume(sCuAnode, Cu, "CuAnode_LV");
  auto lvPEEKBase   = new G4LogicalVolume(sPEEKBase, PEEK, "PEEKBase_LV");
  auto lvPEEKGroove = new G4LogicalVolume(sPEEKGroove, PEEK, "PEEKGroove_LV");
  auto lvTiFrame    = new G4LogicalVolume(sTiFrame, Ti, "TiFrame_LV");
  auto lvBeWin      = new G4LogicalVolume(sBeWin, Be, "BeWindow_LV");

  auto lvTimepix    = new G4LogicalVolume(sTimepix, Si, "TIME_PIX_LV");
  auto lvFoil       = new G4LogicalVolume(sFoil, Al, "Foil_LV");

  // IMPORTANT: this is the SD+field volume
  gas_chamber_Log   = new G4LogicalVolume(sGas, mix, "GAS_CHAMBER_LV");
  gas_sensitive_Log = new G4LogicalVolume(sGasSens, mix, "GAS_SENSITIVE_LV");


  // ---------- VIS ----------
  lvPCBSupport->SetVisAttributes(visWhite);
  lvPCB->SetVisAttributes(visGrey);
  lvVespelBase->SetVisAttributes(visWhite);
  lvVespelGro->SetVisAttributes(visWhite);
  lvCuAnode->SetVisAttributes(visBrown);
  lvPEEKBase->SetVisAttributes(visGreen);
  lvPEEKGroove->SetVisAttributes(visGreen);
  lvTiFrame->SetVisAttributes(visGrey);
  lvBeWin->SetVisAttributes(visGrey);
  gas_chamber_Log->SetVisAttributes(visGold);
  lvTimepix->SetVisAttributes(visGrey);
  lvFoil->SetVisAttributes(visRed);
  gas_sensitive_Log->SetVisAttributes(visSens);

  // ---------- PLACEMENT (single z cursor, +gaps) ----------
  // z0 is the bottom face of PCBSupport.
  G4double z = 0.0;
  const bool overlapCheck = false;

  auto placeAt = [&](G4LogicalVolume* lv, const G4String& name, G4double halfT, G4int copyNo = 0)
  {
    const G4double zc = z + halfT;
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,zc), lv, name, worldLV, false, copyNo, overlapCheck);
    z = zc + halfT + gap;
    return zc;
  };

  // 1) Al support
  const G4double zPCBSupport = placeAt(lvPCBSupport, "PCBSupportAl", PCBSupportAlThick/2., 10);

  // 2) PCB
  const G4double zPCB = placeAt(lvPCB, "PCB", PCB_thick/2., 11);

  // 3) Vespel base
  const G4double zVespelBase = placeAt(lvVespelBase, "VespelBase", Vespel_base_thick/2., 12);

  // 4) Vespel groove ring (same thickness as anode, will sit on top)
  const G4double zVespelGroove = placeAt(lvVespelGro, "VespelGroove", Vespel_groove_thick/2., 13);

  // 5) Cu anode disk: put it at same z as groove (radially disjoint, so no overlap)
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zVespelGroove), lvCuAnode, "CuAnode", worldLV, false, 14, overlapCheck);

  // 6) PEEK base ring (drift length)
  const G4double zPEEKBase = placeAt(lvPEEKBase, "PEEKBase", PEEK_base_thick/2., 15);

  // 7) GAS drift cylinder: place at same z as PEEK base (radially disjoint, so no overlap)
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zPEEKBase), gas_chamber_Log, "GAS_CHAMBER_PV", worldLV, false, 0, overlapCheck);
  // Place the sensitive gas box inside the outer gas cylinder (local coords of GAS_CHAMBER_LV)
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0,0,0),
                    gas_sensitive_Log,
                    "GAS_SENSITIVE_PV",
                    gas_chamber_Log,   // <-- mother is the gas cylinder
                    false,
                    1,
                    overlapCheck);

  // 8) PEEK groove ring
  const G4double zPEEKGroove = placeAt(lvPEEKGroove, "PEEKGroove", PEEK_groove_thick/2., 16);

  // 9) Ti frame disk: same z as groove (radially disjoint)
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zPEEKGroove), lvTiFrame, "TiFrame", worldLV, false, 17, overlapCheck);

  // 10) Be window: place just BELOW Ti frame with a gap (no overlap)
  const G4double zBe = zPEEKGroove - (PEEK_groove_thick/2. + gap + Be_window_thick/2.);
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zBe), lvBeWin, "BeWindow", worldLV, false, 18, overlapCheck);

  // ---------- TIMEPIX + FOIL + PILLARS (tie them to the anode plane) ----------
  // Put Timepix just BELOW the anode plane, then pillars+foil above it.
  const G4double zAnode = zVespelGroove;

  const G4double zTimepix = zAnode - (Cu_Anode_disk_thick/2. + gap + thickness_timepix/2.);
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zTimepix), lvTimepix, "TIMEPIX_PV", worldLV, false, 1, overlapCheck);

  const G4double zFoil = zTimepix + (thickness_timepix/2. + gap + pillarHeight + foilThickness/2.);
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,zFoil), lvFoil, "Al_FOIL_PV", worldLV, false, 56, overlapCheck);

  // Pillars between timepix and foil: center at zTimepix + thickness_timepix/2 + gap + pillarHeight/2
  const G4double zPillar = zTimepix + (thickness_timepix/2. + gap + pillarHeight/2.);

  auto rot = new G4RotationMatrix();
  rot->rotateY(0.0*deg);

  for (int i = 0; i < numPillarsX; ++i) {
    for (int j = 0; j < numPillarsY; ++j) {

      const G4double x = -foilSizeXY/2. + (i + 0.5) * pillarSpacing;
      const G4double y = -foilSizeXY/2. + (j + 0.5) * pillarSpacing;

      auto sPillar = new G4Tubs("Pillar", 0.0, pillarDiameter/2.0, pillarHeight/2.0, 0.0, 360.0*deg);
      G4Transform3D tr(*rot, G4ThreeVector(0,0,0));
      auto sPillarDisp = new G4DisplacedSolid("PillarDisp", sPillar, tr);

      auto lvPillar = new G4LogicalVolume(sPillarDisp, SU8, "Pillar_LV");
      lvPillar->SetVisAttributes(visBlue);

      new G4PVPlacement(nullptr, G4ThreeVector(x,y,zPillar), lvPillar, "Pillar_PV", worldLV, false, 55, overlapCheck);
    }
  }

  // ---------- COMPUTE DETECTOR Z EXTENT ----------
  const G4double zMinDet = 0.0*mm;   // by construction, your first volume starts at z=0
  const G4double zMaxDet = z - gap;  // 'z' is advanced to (top + gap). Remove last gap.
  const G4double detHalfZ = (zMaxDet - zMinDet)/2.0;
  const G4double detZc    = (zMaxDet + zMinDet)/2.0;

  // ---------- GAGG SHIELD SOLIDS ----------
  // Cylindrical wall: hollow tube, open at +Z and -Z (cap handled separately)
  auto sGAGGWall = new G4Tubs("GAGG_Wall_S",
                              shieldRin,   // inner radius
                              shieldRout,  // outer radius
                              detHalfZ,    // half-length along Z
                              0.*deg, 360.*deg);

  // Endcap on -Z side: SOLID disk (no hole)
  // (inner radius = 0 => full disk)
  auto sGAGGCap = new G4Tubs("GAGG_Cap_S",
                            0.,          // inner radius: solid disk
                            shieldRout,  // outer radius
                            shieldThick/2.,  // half-thickness
                            0.*deg, 360.*deg);

  // ---------- GAGG SHIELD LOGICAL VOLUMES ----------
  lvGAGGWall = new G4LogicalVolume(sGAGGWall, GAGG, "GAGG_Wall_LV");
  lvGAGGCap  = new G4LogicalVolume(sGAGGCap,  GAGG, "GAGG_Cap_LV");

  // Make them visible
  auto visGAGG = new G4VisAttributes(G4Colour(0.7, 0.2, 0.8, 0.35)); // purple-ish, semi-transparent
  visGAGG->SetForceSolid(true);
  lvGAGGWall->SetVisAttributes(visGAGG);
  lvGAGGCap->SetVisAttributes(visGAGG);

  // ---------- GAGG SHIELD PLACEMENT ----------
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0,0,detZc),
                    lvGAGGWall,
                    "GAGG_Wall_PV",
                    worldLV,
                    false, 200,
                    overlapCheck);

  // Place cap below z=0 (closed on the -Z side). No +Z cap -> entrance open.
  const G4double zCapC = zMinDet - shieldThick/2.0;
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0,0,zCapC),
                    lvGAGGCap,
                    "GAGG_Cap_PV",
                    worldLV,
                    false, 201,
                    overlapCheck);
  

  // ---------- AL BOX SHELL SOLIDS ----------
  auto sAlOuter = new G4Box("AlBoxOuter_S",
                            alBoxOuterXY/2., alBoxOuterXY/2., alBoxOuterZ/2.);

  auto sAlInner = new G4Box("AlBoxInner_S",
                            (alBoxOuterXY/2. - alBoxThick),
                            (alBoxOuterXY/2. - alBoxThick),
                            (alBoxOuterZ/2.  - alBoxThick) );

  // Hollow shell
  auto sAlShell = new G4SubtractionSolid("AlBoxShell_S", sAlOuter, sAlInner);

  // ---------- +Z ENTRANCE HOLE (circular) ----------
  const G4double holeHalfZ = alBoxThick/2. + 0.1*mm; // tiny margin

  auto sHole = new G4Tubs("AlBoxHole_S",
                          0., alHoleRadius,
                          holeHalfZ,
                          0.*deg, 360.*deg);

  // Place hole cutter at the center of the +Z wall
  const G4double zHoleLocal = +alBoxOuterZ/2. - alBoxThick/2.;

  auto sAlShellWithHole = new G4SubtractionSolid("AlBoxShellWithHole_S",
                                                sAlShell,
                                                sHole,
                                                nullptr,
                                                G4ThreeVector(0,0,zHoleLocal));

  // ---------- AL BOX LOGICAL VOLUME ----------
  auto lvAlBox = new G4LogicalVolume(sAlShellWithHole, Al, "AlBoxShell_LV");

  // visualization
  auto visAlBox = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.15));
  visAlBox->SetForceSolid(true);
  lvAlBox->SetVisAttributes(visAlBox);

  // ---------- PLACEMENT ----------
  new G4PVPlacement(nullptr,
                    G4ThreeVector(0,0,detZc),
                    lvAlBox,
                    "AlBoxShell_PV",
                    worldLV,
                    false, 300,
                    overlapCheck);



  // Print material table (optional)
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return worldPV;
}




/*!
 */
void GPD3D_DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors
  G4String gasCellSDname = "gpd/GasCellSD";
  GPD3D_GasCellSD* gasCellSD = new GPD3D_GasCellSD("gpd/GasCellSD", "GasCellHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(gasCellSD);
  // GPD3D_GasCellSD* GAGG_Wall_SD = new GPD3D_GasCellSD("gpd/GAGG_Wall", "GAGGWallHitsCollection");
  // G4SDManager::GetSDMpointer()->AddNewDetector(GAGG_Wall_SD);
  // GPD3D_GasCellSD* GAGG_Cap_SD = new GPD3D_GasCellSD("gpd/GAGG_Cap", "GAGGCapHitsCollection");
  // G4SDManager::GetSDMpointer()->AddNewDetector(GAGG_Cap_SD);

  // Setting GasCellSD to the gas cell
  SetSensitiveDetector("GAS_SENSITIVE_LV", gasCellSD, true);

  // SetSensitiveDetector("GAGG_Wall_LV", GAGG_Wall_SD, true);
  // SetSensitiveDetector("GAGG_Cap_LV", GAGG_Cap_SD, true);

  // Create the electric field
  GPD3D_DetectorConstruction::ElectricFieldSetup();

}


/*!
 */
void GPD3D_DetectorConstruction::ElectricFieldSetup()
{
  // G4double fieldValue = 2.*kilovolt/cm;
  G4double fieldValue = 1.85*kilovolt/cm;
  m_EMfield = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, fieldValue));

  // Create an equation of motion for this field
  m_equation = new G4EqMagElectricField(m_EMfield);

  G4int nvar = 8;
  m_stepper = new G4ClassicalRK4(m_equation, nvar);

  // Get the local field manager
  G4FieldManager* localFieldManager = new G4FieldManager();

  // Set this field to the local field manager
  localFieldManager->SetDetectorField(m_EMfield);

  m_minStep = 5*um;
  G4MagInt_Driver* intgrDriver = new G4MagInt_Driver(m_minStep, m_stepper,
                                            m_stepper->GetNumberOfVariables());
  m_chordFinder = new G4ChordFinder(intgrDriver);
  localFieldManager->SetChordFinder(m_chordFinder);

  gas_chamber_Log->SetFieldManager(localFieldManager, true);
}


