/**
 * @file
 * @brief Default module that is always loaded by Core
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "MrccIO.h"
#include "Utils/IO/NativeFilenames.h"
#include <Utils/CalculatorBasics/CalculationRoutines.h>
#include <Utils/Constants.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/ExternalQC/SettingsNames.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/Regex.h>
#include <Utils/Strings.h>
#include <Utils/UniversalSettings/SettingsNames.h>
#include <boost/algorithm/string.hpp>
#include <boost/process.hpp>
#include <cmath>
#include <iomanip>
#include <regex>
#include <vector>

namespace Scine {
namespace Utils {
namespace ExternalQC {

MrccIO::MrccIO(const MrccFiles& files, const Settings& settings, const std::string& methodFamily)
  : files_(files), settings_(settings), mrccMethod_(getMrccMethod(settings, methodFamily)) {
}

void MrccIO::writeInput(const AtomCollection& atoms) {
  std::ofstream outStream;
  outStream.open(files_.input);
  this->addAllowedResources(outStream);
  this->addChargeAndMultiplicity(outStream);
  this->addMethodDefinition(outStream);
  this->addSCFKeywords(outStream);
  this->addBasisSetKeyword(outStream);
  this->addSolvationKeywords(outStream);
  this->addSCFTypeKeyword(outStream);
  // The coordinate definition should remain the last part of the input because we just append the xyz-file.
  this->addCoordinateDefinition(atoms, outStream);
  outStream.close();
}

void MrccIO::addAllowedResources(std::ofstream& outStream) {
  outStream << "mem=" << settings_.getInt(Utils::SettingsNames::externalProgramMemory) << "mb" << std::endl;
}

void MrccIO::addChargeAndMultiplicity(std::ofstream& outStream) {
  outStream << "mult=" << settings_.getInt(Utils::SettingsNames::spinMultiplicity) << std::endl;
  outStream << "charge=" << settings_.getInt(Utils::SettingsNames::molecularCharge) << std::endl;
}

void MrccIO::addMethodDefinition(std::ofstream& outStream) {
  this->addCalcKeyword(outStream);
  if (this->isLocalCorrelation()) {
    this->addLocalCorrelationKeywords(outStream);
  }
}

void MrccIO::addLocalCorrelationKeywords(std::ofstream& outStream) {
  outStream << "lcorthr=" << this->getLNOThresholdsFromMethod() << std::endl;
  outStream << "core=frozen" << std::endl;
  outStream << "ccsalg=dfdirect" << std::endl;
  outStream << "ccprog=ccsd" << std::endl;
}

void MrccIO::addCalcKeyword(std::ofstream& outStream) {
  switch (mrccMethod_) {
    case MrccMethod::HF: {
      outStream << "calc=hf" << std::endl;
      return;
    }
    case MrccMethod::DFT: {
      outStream << "calc=" << this->functionalInMrccFormat() << std::endl;
      return;
    }
    case MrccMethod::LNO_MP2: {
      outStream << "calc=lno-mp2" << std::endl;
      return;
    }
    case MrccMethod::LNO_CCSD: {
      outStream << "calc=lno-ccsd" << std::endl;
      outStream << "localcc=on" << std::endl;
      return;
    }
    case MrccMethod::LNO_CCSD_T: {
      outStream << "calc=lno-ccsd(t)" << std::endl;
      outStream << "localcc=on" << std::endl;
      return;
    }
  }
  throw std::logic_error("MRCC method not handled in switch.");
}

void MrccIO::addSCFTypeKeyword(std::ofstream& outStream) {
  const auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings_.getString(Utils::SettingsNames::spinMode));
  switch (spinMode) {
    case SpinMode::None:
    case SpinMode::Any:
      return;
    case SpinMode::Restricted: {
      outStream << "scftype=RHF" << std::endl;
      return;
    }
    case SpinMode::Unrestricted: {
      outStream << "scftype=UHF" << std::endl;
      return;
    }
    case SpinMode::RestrictedOpenShell: {
      outStream << "scftype=ROHF" << std::endl;
      return;
    }
  }
  throw std::runtime_error("Spin mode not handled in MRCC calculator.");
}

void MrccIO::addSCFKeywords(std::ofstream& outStream) {
  if (settings_.getBool(Utils::SettingsNames::scfDamping))
    outStream << "scfdamp=" << settings_.getDouble(SettingsNames::scfDampingValue) << std::endl;
  outStream << "scflshift=" << settings_.getDouble(SettingsNames::scfOrbitalShift) << std::endl;
  const double scfTol = settings_.getDouble(Utils::SettingsNames::selfConsistenceCriterion);
  outStream << "scftol=" << int(round(-log10(scfTol))) << std::endl;
  outStream << "scfmaxit=100" << std::endl;
}

void MrccIO::addSolvationKeywords(std::ofstream& outStream) {
  const std::string solvation = settings_.getString(Utils::SettingsNames::solvation);
  const std::string solvent = settings_.getString(Utils::SettingsNames::solvent);
  if (solvation == "iefpcm") {
    outStream << "pcm=" << solvent << std::endl;
  }
}

void MrccIO::addBasisSetKeyword(std::ofstream& outStream) {
  const std::string basis = settings_.getString(Utils::SettingsNames::basisSet);
  outStream << "basis=" << basis << std::endl;
}

void MrccIO::addCoordinateDefinition(const AtomCollection& atoms, std::ofstream& outStream) {
  outStream << "geom=xyz" << std::endl;
  outStream << atoms.size() << "\n" << std::endl;
  for (const auto& atom : atoms) {
    std::string elementName = ElementInfo::symbol(atom.getElementType());
    outStream << std::setw(4) << std::left;
    outStream << elementName << atom.getPosition() * Constants::angstrom_per_bohr << "\n";
  }
  outStream << std::right << std::setw(0);
  outStream.flush();
}

std::string MrccIO::getLNOThresholdsFromMethod() {
  std::string method = settings_.getString(Utils::SettingsNames::method);
  boost::algorithm::to_lower(method);
  std::vector<std::string> lnoOptions = {"vloose", "loose", "normal", "tight", "vtight", "vvtight"};
  for (const auto& option : lnoOptions) {
    if (method.find(option) != std::string::npos)
      return option;
  }
  this->getLog().warning << "No LNO threshold definition detected for the local correlation calculation with MRCC."
                         << " The calculation will be performed with 'normal' settings."
                         << " Input example: tight-lno-ccsd(t)" << Core::Log::nl;
  return "normal";
}

bool MrccIO::isLocalCorrelation() {
  return mrccMethod_ == MrccMethod::LNO_MP2 || mrccMethod_ == MrccMethod::LNO_CCSD || mrccMethod_ == MrccMethod::LNO_CCSD_T;
}

std::string MrccIO::readOutput() {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  if (!boost::filesystem::exists(this->files_.output)) {
    throw std::runtime_error("File " + this->files_.output + " not found.");
  }
  fin.open(this->files_.output);
  std::string out = std::string(std::istreambuf_iterator<char>{fin}, {});
  fin.close();
  ensureSuccessFullCalculation(out);
  return out;
}

std::string MrccIO::getEnergyString() {
  switch (mrccMethod_) {
    case MrccMethod::HF:
      return "FINAL HARTREE-FOCK ENERGY:";
    case MrccMethod::DFT:
      return this->functionalInMrccFormat() + " energy \\[au\\]:";
    case MrccMethod::LNO_MP2:
      return R"(DF-MP2 energy \[au\]:)";
    case MrccMethod::LNO_CCSD:
      return R"(Total LNO-CCSD energy with MP2 corrections \[au\]:)";
    case MrccMethod::LNO_CCSD_T:
      return R"(Total LNO-CCSD\(T\) energy with MP2 corrections \[au\]:)";
  }
  throw std::logic_error("MRCC method option not handled in switch.");
}

double MrccIO::getEnergy(const std::string& output) {
  std::string regexString = this->getEnergyString() + " + " + Regex::capturingFloatingPointNumber();
  std::regex regex(regexString);
  std::sregex_iterator iter(output.begin(), output.end(), regex);
  std::sregex_iterator end;
  double energy = 0.0;
  bool found = false;
  while (iter != end) {
    energy = std::stod((*iter)[1]);
    found = true;
    ++iter;
  }
  if (!found) {
    throw OutputFileParsingError("Energy could not be read from MRCC output.");
  }
  return energy;
}

void MrccIO::ensureSuccessFullCalculation(const std::string& out) {
  const std::string scfFailure = "THE SCF ITERATION HAS NOT CONVERGED";
  const std::string normalTermination = "Normal termination of mrcc";
  if (out.find(scfFailure) != std::string::npos) {
    throw ScfNotConvergedError("MRCC could not converge the SCF.");
  }
  if (out.find(normalTermination) == std::string::npos) {
    throw OutputFileParsingError("MRCC did not terminate normally.");
  }
}

MrccMethod MrccIO::getMrccMethod(const Settings& settings, const std::string& methodFamily) {
  if (caseInsensitiveEqual(methodFamily, "hf")) {
    return MrccMethod::HF;
  }
  if (caseInsensitiveEqual(methodFamily, "dft")) {
    return MrccMethod::DFT;
  }
  if (caseInsensitiveEqual(methodFamily, "mp2")) {
    return MrccMethod::LNO_MP2;
  }
  if (caseInsensitiveEqual(methodFamily, "cc")) {
    std::string method = settings.getString(Utils::SettingsNames::method);
    boost::algorithm::to_lower(method);
    if (method.find("ccsd(t)") != std::string::npos) {
      return MrccMethod::LNO_CCSD_T;
    }
    if (method.find("ccsd") != std::string::npos) {
      return MrccMethod::LNO_CCSD;
    }
    throw std::runtime_error("Electronic structure method " + method + " combined with the method family " +
                             methodFamily + "  is unknown to MRCC.");
  }
  throw std::runtime_error("Unknown electronic structure method family to MRCC " + methodFamily);
}

std::string MrccIO::functionalInMrccFormat() {
  const std::string functional = settings_.getString(Utils::SettingsNames::method);
  auto methodInput = Scine::Utils::CalculationRoutines::splitIntoMethodAndDispersion(functional);
  boost::to_upper(methodInput.first);
  if (methodInput.second.empty()) {
    return methodInput.first;
  }
  if (caseInsensitiveEqual("D3BJ", methodInput.second)) {
    return methodInput.first + "-D3";
  }
  throw std::runtime_error("The SCINE-MRCC interface supports only D3BJ as dispersion correction.");
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
