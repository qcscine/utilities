/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "TurbomoleOrbitalSteerer.h"
#include "TurbomoleOrbitalFileReader.h"
#include "TurbomoleOrbitalFileWriter.h"
#include "Utils/IO/FilesystemHelpers.h"
#include <Core/Log.h>
#include <Utils/DataStructures/MolecularOrbitals.h>
#include <Utils/IO/NativeFilenames.h>
#include <Utils/Scf/OrbitalPerturbation/RandomOrbitalMixer.h>
#include <regex>
#include <string>

namespace Scine {
namespace Utils {
namespace ExternalQC {
namespace TurbomoleOrbitalPerturbation {

TurbomoleOrbitalSteerer::TurbomoleOrbitalSteerer(std::string& calculationDirectory)
  : calculationDirectory_(calculationDirectory) {
}

void TurbomoleOrbitalSteerer::steerOrbitals() {
  setCorrectTurbomoleFileNames(files_, calculationDirectory_);
  std::pair<int, int> nElec = extractNumberOfAlphaAndBetaElectrons();
  int nOrbitals = extractNumerOfOrbitals();

  TurbomoleOrbitalFileReader ta(files_.alphaFile, nOrbitals);
  TurbomoleOrbitalFileReader tb(files_.betaFile, nOrbitals);

  FilesystemHelpers::copyFile(files_.alphaFile, files_.alphaBakFile);
  FilesystemHelpers::copyFile(files_.betaFile, files_.betaBakFile);

  MolecularOrbitals mo =
      MolecularOrbitals::createFromUnrestrictedCoefficients(ta.getCoefficientMatrix(), tb.getCoefficientMatrix());
  mixOrbitals(mo, nElec.first, nElec.second);
  TurbomoleOrbitalFileWriter wa(mo.alphaMatrix(), ta.getMetaInformation());
  TurbomoleOrbitalFileWriter wb(mo.betaMatrix(), tb.getMetaInformation());
  wa.writeToFile(files_.alphaFile);
  wb.writeToFile(files_.betaFile);
}

void TurbomoleOrbitalSteerer::mixOrbitals(MolecularOrbitals& mo, unsigned nAlpha, unsigned nBeta) {
  OrbitalPerturbation::RandomOrbitalMixer mixer(mo, nAlpha, nBeta);
  Core::Log log;
  mixer.setMaximalMixAngle(1.57);
  mixer.setNumberMixes(10);
  mixer.considerOnlyOrbitalsCloseToFrontierOrbitals(15);
  mixer.mix(log);
}

std::pair<int, int> TurbomoleOrbitalSteerer::extractNumberOfAlphaAndBetaElectrons() {
  int nAlpha = 0;
  int nBeta = 0;
  std::ifstream in;
  in.open(NativeFilenames::combinePathSegments(calculationDirectory_, "control"));
  std::string line;
  std::regex r(R"([a]+ +([0-9]+)(-)([0-9]+))");
  std::smatch m1, m2;
  while (std::getline(in, line)) {
    if (line.find("$alpha shells") != std::string::npos) {
      std::string nextLine;
      std::getline(in, nextLine);
      if (std::regex_search(nextLine, m1, r))
        nAlpha = std::stoi(m1[3]);
      else
        throw std::runtime_error("Could not parse number of alpha electrons.");
    }
    else if (line.find("$beta shells") != std::string::npos) {
      std::string nextLine;
      std::getline(in, nextLine);
      if (std::regex_search(nextLine, m2, r))
        nBeta = std::stoi(m2[3]);
      else
        throw std::runtime_error("Could not parse number of beta electrons.");
    }
  }

  in.close();

  return std::make_pair(nAlpha, nBeta);
}

int TurbomoleOrbitalSteerer::extractNumerOfOrbitals() {
  int nOrbitals = 0;
  std::ifstream in;
  in.open(files_.controlFile);
  std::string line;
  std::regex r(R"(nbf\(AO\)=([0-9]+))");
  std::smatch m1;
  while (std::getline(in, line)) {
    if (line.find("$rundimensions") != std::string::npos) {
      std::string nextLine;
      std::getline(in, nextLine);
      std::getline(in, nextLine);
      std::getline(in, nextLine);
      if (std::regex_search(nextLine, m1, r))
        nOrbitals = std::stoi(m1[1]);
      else
        throw std::runtime_error("Could not parse number of orbitals.");
    }
  }

  in.close();

  return nOrbitals;
}

} // namespace TurbomoleOrbitalPerturbation
} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
