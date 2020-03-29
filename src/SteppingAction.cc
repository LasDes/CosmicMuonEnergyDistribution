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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction *eventAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *step)
{
  // get detector geometry
  const DetectorConstruction *detectorConstruction = static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  // detector volume
  auto detVolume = detectorConstruction->GetScoringVolume();

  //G4cout << detVolume->GetLogicalVolume()->GetMass() / CLHEP::kg << G4endl;

  // what to look out for?
  // what volume are we in?
  //G4LogicalVolume *volume = step->GetTrack()->GetVolume()->GetLogicalVolume();
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // what kind of particle are we talking about?
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

  // open file and write to file
  std::ofstream data;

  // track length in detector for muons
  //if (volume == detectorConstruction->GetScoringVolume() && particleName == "mu-" && trackID == 1)
  if (volume == detVolume)
  {
    G4double stepLength = step->GetStepLength() / CLHEP::mm;
    G4double edep = step->GetTotalEnergyDeposit() / CLHEP::keV;
    fEventAction->AddTrackLength(stepLength);
    fEventAction->AddMuonEnergy(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
