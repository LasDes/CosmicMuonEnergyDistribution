//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fScoringVolume(0),
      fCase(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    // Get nist material manager
    G4NistManager *nist = G4NistManager::Instance();

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4double world_sizeXY = 10 * m;
    G4double world_sizeZ = 10 * m;
    G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box *solidWorld = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);

    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

    //
    // detector zone
    //

    // Al case
    // sizes
    G4double AlCaseRadius = 3.8 * cm;
    G4double AlCaseHeight = 50. / 2. * cm;

    // material --> Al
    G4Material *caseMaterial = nist->FindOrBuildMaterial("G4_Al");

    // solid
    G4Tubs *caseSolid = new G4Tubs("AlCase", 0, AlCaseRadius, AlCaseHeight, 0, 2 * M_PI);

    // inner space
    // sizes
    G4double innerRadius = 3.6 * cm;
    G4double innerHeight = 49.6 / 2. * cm;

    // material P10(?) (10% CO2, 90% Ar)
    G4int num_elements;
    G4double density = 1.8033 * kg / m3;
    G4Material *elAr = nist->FindOrBuildMaterial("G4_Ar");
    G4Material *matCO2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    G4Material *innerMaterial = new G4Material("innerMaterial", density, num_elements = 2, kStateGas, 300. * kelvin, 1 * atmosphere);
    innerMaterial->AddMaterial(matCO2, .1);
    innerMaterial->AddMaterial(elAr, .9);

    // solid
    G4Tubs *innerSolid = new G4Tubs("innerSpace", 0, innerRadius, innerHeight, 0, 2 * M_PI);

    // cut smaller cylinder from bigger
    G4SubtractionSolid *hollowCase = new G4SubtractionSolid("hollowCase", caseSolid, innerSolid);

    // logical volume (case)
    G4LogicalVolume *caseLogical = new G4LogicalVolume(hollowCase, caseMaterial, "AlCase");
    // physical volume (case)
    new G4PVPlacement(0, G4ThreeVector(), caseLogical, "AlCase", logicWorld, false, 0, checkOverlaps);

    // logical volume (inner)
    G4LogicalVolume *innerLogical = new G4LogicalVolume(innerSolid, innerMaterial, "innerSpace");
    // physical volume (inner)
    G4VPhysicalVolume *innerPhysical = new G4PVPlacement(0, G4ThreeVector(), innerLogical, "innerSpace", logicWorld, false, 0, checkOverlaps);

    //make inner space the scoring volume
    fScoringVolume = innerPhysical;

    //
    //always return the physical World
    //
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
