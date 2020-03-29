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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

// std random number generation
#include <random>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // set particle parameters
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "mu-");
  fParticleGun->SetParticleDefinition(particle);

  // particle initial position
  G4double z = 25. * cm;
  G4ThreeVector pos(10. * cm, 0., z);
  fParticleGun->SetParticlePosition(pos);

  // paticle energy
  fParticleGun->SetParticleEnergy(4. * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// inverse function approximation --> bisection
template <typename T, typename TT, typename TTT, typename F>
T Inverse_Func_Approximation(T x, TT low, TTT high, F func)
{
  T eps = 1e-15;
  while ((high - low) > eps)
  {
    T mid = (low + high) / 2;
    if (func(mid) == x)
    {
      return mid;
    }
    else if (func(mid) < x)
    {
      low = mid;
    }
    else
    {
      high = mid;
    }
  }

  return (low + high) / 2;
}

// muons' angular cummulative distribution function --> lambda
auto Muonic_Angular_CDF = [](G4double x) { return (2 * x + std::sin(2 * x) - pi) / pi; };
// square function --> lambda
auto sq = [](G4double x) { return x * x; };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  // this function is called at the begining of each event
  // useful parameters
  // must be updated if detector size changes in DetectorConstruction... poor hard code :(
  // height of full detector
  G4double Height = 50. * cm;
  // radius of full detector
  G4double Radius = 3.8 * cm;

  // position vector
  G4ThreeVector pos = fParticleGun->GetParticlePosition();
  // distance from origin (2D polar coordianates)
  G4double r = std::sqrt(sq(pos.getX()) + sq(pos.getY()));

  // initialize zenith and azimuth angle
  G4double theta = pi;
  G4double phi = 0.;

  // seperate cases for distances
  // in the detector radius
  if (r <= Radius)
  {
    // theta zenith angle for momentum direction is arbitrary from [pi / 2, pi] with cos^2(theta) distribution
    theta = Inverse_Func_Approximation(G4UniformRand(), halfpi, pi, Muonic_Angular_CDF);
    
    // phi azimuth angle for momentum direction is arbitrary from [0, 2pi] with uniform distibution
    phi = G4UniformRand() * twopi;
  }
  // out of the detector radius
  else if (r > Radius)
  {
    G4double phiMin = 0., phiMax = 0.;

    // importat angles
    G4double thetaB = std::atan((r - Radius) / Height) + halfpi, thetaMax = std::atan(std::sqrt(sq(r) + sq(Radius)) / Height) + halfpi;

    // theta zenith angle for momentum direction is arbitrary from [pi / 2, pi] with cos^2(theta) distribution
    while (theta > thetaMax)
    {
      theta = Inverse_Func_Approximation(G4UniformRand(), halfpi, pi, Muonic_Angular_CDF);
    }

    // further checks for theta to determine boundaries for the possible azimuth angle inteval
    if (theta <= thetaB)
    {
      G4double delta0 = std::asin(Radius / r);
      phiMin = pi - delta0, phiMax = pi + delta0;
    }
    else if ((theta > thetaB) && (theta <= thetaMax))
    {
      // new variable xi = theta - pi / 2
      G4double tanXi = std::tan(theta - halfpi);
      G4double delta = std::asin(Radius / Height / tanXi * std::sqrt(1 - sq(sq(r) + sq(Radius) - sq(Height * tanXi)) / sq(2 * r * Radius)));
      phiMin = pi - delta, phiMax = pi + delta;
    }

    // random number generation via std
    std::random_device rd{};
    std::mt19937 gen(rd());
    std::uniform_real_distribution<G4double> distributionPhi(phiMin, phiMax);
    phi = distributionPhi(gen);
  }

  // test distance
  G4cout << r / cm << " cm" << G4endl;

  // setting momentum direction vector
  G4double sinTheta = std::sin(theta);
  G4double vx = sinTheta * std::cos(phi), vy = sinTheta * std::sin(phi), vz = std::cos(theta);
  G4ThreeVector momentumDir(vx, vy, vz);

  fParticleGun->SetParticleMomentumDirection(momentumDir);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
