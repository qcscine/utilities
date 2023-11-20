/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Cp2kInputFileCreator.h"
#include "Cp2kCalculatorSettings.h"
#include <Utils/CalculatorBasics.h>
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <algorithm>
#include <fstream>
#include <iterator>

namespace Scine {
namespace Utils {
namespace ExternalQC {

Cp2kInputFileCreator::Cp2kInputFileCreator(const AtomCollection& atoms, const Settings& settings,
                                           const PropertyList& requiredProperties, bool isDft)
  : _atoms(atoms), _settings(settings), _requiredProperties(requiredProperties), _isDft(isDft) {
}

void Cp2kInputFileCreator::createInputFile(const std::string& filename, const std::string& projectName) const {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, projectName);
  fout.close();
  CalculationRoutines::checkValidityOfChargeAndMultiplicity(_settings.getInt(Utils::SettingsNames::molecularCharge),
                                                            _settings.getInt(Utils::SettingsNames::spinMultiplicity), _atoms);
}

void Cp2kInputFileCreator::createInputFile(std::ostream& out, const std::string& projectName) const {
  printGlobal(out, projectName);
  printForceEval(out);
}

void Cp2kInputFileCreator::printGlobal(std::ostream& out, const std::string& projectName) const {
  std::string runType = (_requiredProperties.containsSubSet(Property::Hessian)) ? "VIBRATIONAL_ANALYSIS" : "ENERGY_FORCE";
  std::string printLevel = "MEDIUM";
  out << "&GLOBAL" << std::endl;
  out << "\tPROJECT " << projectName << std::endl;
  out << "\tRUN_TYPE " << runType << std::endl;
  out << "\tPRINT_LEVEL " << printLevel << std::endl;
  out << "\tEXTENDED_FFT_LENGTHS" << std::endl; // necessary for slab calculations
  out << "\t&PRINT" << std::endl;
  out << "\t\tPHYSCON FALSE" << std::endl; // no physical constants
  out << "\t&END PRINT" << std::endl;
  out << "&END GLOBAL" << std::endl << std::endl;
}

void Cp2kInputFileCreator::printForceEval(std::ostream& out) const {
  out << "&FORCE_EVAL" << std::endl;
  out << "\t&PRINT" << std::endl;
  out << "\t\t&FORCES ON" << std::endl;
  out << "\t\t\tNDIGITS 18" << std::endl;
  out << "\t\t&END FORCES" << std::endl;
  if (_requiredProperties.containsSubSet(Property::StressTensor)) {
    out << "\t\t&STRESS_TENSOR ON" << std::endl;
    out << "\t\t&END STRESS_TENSOR" << std::endl;
  }
  out << "\t&END PRINT" << std::endl;
  out << "\tMETHOD QUICKSTEP" << std::endl;
  if (_requiredProperties.containsSubSet(Property::StressTensor)) {
    out << "\tSTRESS_TENSOR ANALYTICAL" << std::endl;
  }
  printSubsys(out);
  // add method variations here if something other than DFT has been implemented
  printDftInput(out);
  out << "&END FORCE_EVAL" << std::endl;
}

void Cp2kInputFileCreator::printSubsys(std::ostream& out) const {
  out << "\t&SUBSYS" << std::endl;
  out << "\t\t&PRINT" << std::endl;
  out << "\t\t\t&SYMMETRY" << std::endl;
  out << "\t\t\t\tALL" << std::endl;
  if (_settings.getString(Utils::SettingsNames::periodicBoundaries) == SettingsNames::moleculePeriodicBoundaries) {
    out << "\t\t\t\tMOLECULE" << std::endl;
  }
  out << "\t\t\t&END SYMMETRY" << std::endl;
  out << "\t\t&END PRINT" << std::endl;
  printCell(out);
  printCoords(out);
  printBasis(out);
  out << "\t&END SUBSYS" << std::endl;
}

void Cp2kInputFileCreator::printCell(std::ostream& out) const {
  out << "\t\t&CELL" << std::endl;
  auto pbc = PeriodicBoundaries(_settings.getString(Utils::SettingsNames::periodicBoundaries));
  out << "\t\t\tABC " << std::to_string(pbc.getA().norm() * Constants::angstrom_per_bohr) << " "
      << std::to_string(pbc.getB().norm() * Constants::angstrom_per_bohr) << " "
      << std::to_string(pbc.getC().norm() * Constants::angstrom_per_bohr) << std::endl;
  out << "\t\t\tALPHA_BETA_GAMMA " << std::to_string(pbc.getAlpha()) << " " << std::to_string(pbc.getBeta()) << " "
      << std::to_string(pbc.getGamma()) << std::endl;
  out << "\t\t\tPERIODIC XYZ" << std::endl;
  out << "\t\t&END CELL" << std::endl;
}

void Cp2kInputFileCreator::printCoords(std::ostream& out) const {
  out << "\t\t&COORD" << std::endl;
  for (const auto& atom : _atoms) {
    out << "\t\t";
    MolecularTrajectoryIO::writeXYZLine(out, atom.getElementType(), atom.getPosition());
  }
  out << "\t\t&END COORD" << std::endl;
  out << "\t\t&TOPOLOGY\n\t\t\t&CENTER_COORDINATES\n\t\t\t&END\n\t\t&END TOPOLOGY" << std::endl; // centering coords
}

void Cp2kInputFileCreator::printBasis(std::ostream& out) const {
  // get unique elements
  auto uniqueElements = _atoms.getElements();
  std::sort(uniqueElements.begin(), uniqueElements.end());
  auto it = std::unique(uniqueElements.begin(), uniqueElements.end());
  uniqueElements.resize(std::distance(uniqueElements.begin(), it));
  if (!_isDft) {
    // write different elements
    for (const auto& element : uniqueElements) {
      out << "\t\t&KIND " << ElementInfo::symbol(element) << std::endl;
      out << "\t\t\tELEMENT " << ElementInfo::symbol(element) << std::endl;
      out << "\t\t&END KIND" << std::endl;
    }
    return; // we do not need basis sets for semi empirical methods
  }
  auto basisEntry = _settings.getString(Utils::SettingsNames::basisSet);
  std::for_each(basisEntry.begin(), basisEntry.end(), [](char& c) { c = ::toupper(c); });
  if (basisEntry.empty() || basisEntry == "none") {
    throw std::logic_error("You need to specify a basis set if you want to do DFT calculations");
  }
  // sanity checks
  if (basisEntry.find("MOLOPT") == std::string::npos) {
    throw std::logic_error("Currently we only support MOLOPT basis sets, you demanded " + basisEntry);
  }
  if (basisEntry.find("GTH") == std::string::npos) {
    throw std::logic_error("Currently we only support MOLOPT basis sets with GTH pseudopotentials, you demanded " + basisEntry);
  }
  if (std::count(basisEntry.begin(), basisEntry.end(), '-') != 2) {
    throw std::logic_error("Invalid basis set entry, please specify like this 'DZVP-MOLOPT-GTH', you specified " + basisEntry);
  }
  std::string zetaInfo = basisEntry.substr(0, basisEntry.find('-'));
  if (std::find(_availableZetaStrings.begin(), _availableZetaStrings.end(), zetaInfo) == _availableZetaStrings.end()) {
    throw std::logic_error("The provided basis set " + basisEntry + " is not supported");
  }
  // check valid basis
  if (zetaInfo.find('T') != std::string::npos) {
    // triple zetas not available for SR basis sets
    for (const auto& element : uniqueElements) {
      if (std::find(_highQualityBasisElements.begin(), _highQualityBasisElements.end(), element) ==
          _highQualityBasisElements.end()) {
        throw std::logic_error("The provided basis set " + basisEntry + " is not supported for element " +
                               ElementInfo::symbol(element));
      }
    }
  }
  // get method for pseudopotential
  // change here if more pseudopotentials have to be added
  auto method = CalculationRoutines::splitIntoMethodAndDispersion(_settings.getString(Utils::SettingsNames::method)).first;
  std::for_each(method.begin(), method.end(), [](char& c) { c = ::toupper(c); });
  if (method.find("PBE") != std::string::npos) {
    // change to plain PBE in case we have REVPBE or PBESOL
    // because no separate pseudopotential exists for them by default
    method = "PBE";
  }
  // write different elements
  for (const auto& element : uniqueElements) {
    std::string actualBasis = zetaInfo + "-MOLOPT";
    if (std::find(_highQualityBasisElements.begin(), _highQualityBasisElements.end(), element) ==
        _highQualityBasisElements.end()) {
      actualBasis += "-SR";
    }
    actualBasis += "-GTH";

    out << "\t\t&KIND " << ElementInfo::symbol(element) << std::endl;
    out << "\t\t\tELEMENT " << ElementInfo::symbol(element) << std::endl;
    out << "\t\t\tBASIS_SET " << actualBasis << std::endl;
    out << "\t\t\tPOTENTIAL GTH-" + method << std::endl;
    out << "\t\t&END KIND" << std::endl;
  }
}

void Cp2kInputFileCreator::printDftInput(std::ostream& out) const {
  // Semiempirics specifications are within DFT section
  out << "\t&DFT" << std::endl;
  printElectronicStructureBasics(out);
  if (_isDft) {
    printFunctional(out);
  }
  else {
    printSemiempiricalMethod(out);
  }
  printScfInput(out);
  printPoissonSolver(out);
  printGridInput(out);
  printMatrixPrint(out);
  out << "\t&END DFT" << std::endl;
}

void Cp2kInputFileCreator::printElectronicStructureBasics(std::ostream& out) const {
  out << "\t\tCHARGE " << _settings.getInt(Utils::SettingsNames::molecularCharge) << std::endl;
  int multiplicity = _settings.getInt(Utils::SettingsNames::spinMultiplicity);
  out << "\t\tMULTIPLICITY " << multiplicity << std::endl;
  auto spinMode = SpinModeInterpreter::getSpinModeFromString(_settings.getString(Utils::SettingsNames::spinMode));
  out << "\t\t" << determineCp2kSpinMode(spinMode, multiplicity) << std::endl;
}

void Cp2kInputFileCreator::printFunctional(std::ostream& out) const {
  out << "\t\tBASIS_SET_FILE_NAME BASIS_MOLOPT" << std::endl; // TODO modify if more basis sets are implemented
  out << "\t\t&XC" << std::endl;
  auto methodInputs = CalculationRoutines::splitIntoMethodAndDispersion(_settings.getString(Utils::SettingsNames::method));
  auto method = methodInputs.first;
  std::for_each(method.begin(), method.end(), [](char& c) { c = ::toupper(c); });
  out << "\t\t\t&XC_FUNCTIONAL ";
  if (method == "REVPBE" || method == "PBESOL") {
    out << "\n\t\t\t\t&PBE" << std::endl;
    out << "\t\t\t\t\tPARAMETRIZATION " << method << std::endl;
    out << "\t\t\t\t&END PBE" << std::endl;
  }
  else {
    out << method << std::endl;
  }
  out << "\t\t\t&END XC_FUNCTIONAL" << std::endl;
  printDispersionCorrection(out, methodInputs);
  out << "\t\t&END XC" << std::endl;
  if (_settings.getBool(SettingsNames::dipoleCorrection)) {
    out << "\t\tSURFACE_DIPOLE_CORRECTION" << std::endl;
  }
}

void Cp2kInputFileCreator::printDispersionCorrection(std::ostream& out, std::pair<std::string, std::string> methodInputs) const {
  auto vdwFunc = methodInputs.second;
  if (!vdwFunc.empty()) {
    // check if valid setting
    if (_dispersionSettingMap.find(vdwFunc) == _dispersionSettingMap.end()) {
      std::string availables = "[";
      for (auto const& pair : _dispersionSettingMap) {
        availables += pair.first + ", ";
      }
      availables = availables.substr(0, availables.size() - 2); // - 2 because we remove last ", "
      availables += "]";
      throw std::logic_error("Invalid dispersion correction setting '" + vdwFunc + "'. Available inputs are:\n" + availables);
    }
    std::string cp2kDispersionInput = _dispersionSettingMap.at(vdwFunc);
    out << "\t\t\t&VDW_POTENTIAL" << std::endl;
    if (cp2kDispersionInput.find("DFTD") != std::string::npos) {
      out << "\t\t\t\tPOTENTIAL_TYPE PAIR_POTENTIAL" << std::endl;
      out << "\t\t\t\t&PAIR_POTENTIAL" << std::endl;
      out << "\t\t\t\t\tTYPE " << cp2kDispersionInput << std::endl;
      out << "\t\t\t\t\tREFERENCE_FUNCTIONAL " << methodInputs.first << std::endl;
      out << "\t\t\t\t\tPARAMETER_FILE_NAME dftd3.dat" << std::endl;
      out << "\t\t\t\t&END PAIR_POTENTIAL" << std::endl;
    }
    else {
      out << "\t\t\t\tPOTENTIAL_TYPE NON_LOCAL" << std::endl;
      out << "\t\t\t\t&NON_LOCAL" << std::endl;
      out << "\t\t\t\t\tTYPE " << cp2kDispersionInput << std::endl;
      out << "\t\t\t\t\tKERNEL_FILE_NAME vdW_kernel_table.dat" << std::endl;
      out << "\t\t\t\t&END NON_LOCAL" << std::endl;
    }
    out << "\t\t\t&END VDW_POTENTIAL" << std::endl;
  }
}

void Cp2kInputFileCreator::printSemiempiricalMethod(std::ostream& out) const {
  auto method = _settings.getString(Utils::SettingsNames::method);
  std::for_each(method.begin(), method.end(), [](char& c) { c = ::toupper(c); });
  if (method != "GFN1") {
    throw std::logic_error("Currently GFN1 is the only supported semiempirical method");
  }
  out << "\t\t&QS" << std::endl;
  out << "\t\t\tMETHOD xTB" << std::endl;
  out << "\t\t\t&XTB" << std::endl;
  out << "\t\t\t\tDO_EWALD T" << std::endl;
  out << "\t\t\t\tCHECK_ATOMIC_CHARGES False" << std::endl; // SCF guess can be bad with xTB causing error
  out << "\t\t\t\t&PARAMETER" << std::endl;
  out << "\t\t\t\t\tDISPERSION_PARAMETER_FILE dftd3.dat" << std::endl;
  out << "\t\t\t\t&END PARAMETER" << std::endl;
  out << "\t\t\t&END XTB" << std::endl;
  out << "\t\t&END QS" << std::endl;
}

std::string Cp2kInputFileCreator::determineCp2kSpinMode(SpinMode spinMode, int multiplicity) const {
  if (spinMode == SpinMode::None) {
    throw std::logic_error("Specified spin mode to 'none'. This is not possible with DFT.");
  }
  if (spinMode == SpinMode::Restricted) {
    if (multiplicity != 1) {
      throw std::logic_error("Specified restricted spin for multiplicity larger than 1.");
    }
    return "!restricted"; // does not exist as option for cp2k, which enforces it implicitly,
                          // this simply adds a comment to the input file
  }
  if (spinMode == SpinMode::Any) {
    return (multiplicity == 1) ? "!restricted" : "UKS";
  }
  if (spinMode == SpinMode::Unrestricted) {
    return "UKS";
  }
  if (spinMode == SpinMode::RestrictedOpenShell) {
    return "ROKS";
  }
  throw std::logic_error("Specified unknown spin mode " + SpinModeInterpreter::getStringFromSpinMode(spinMode) +
                         " in _settings."); // this should have been handled by _settings
}

void Cp2kInputFileCreator::printScfInput(std::ostream& out) const {
  // E accuracy
  out << "\t\t&SCF" << std::endl;
  // single keywords
  out << "\t\t\tSCF_GUESS " << _settings.getString(SettingsNames::scfGuess) << std::endl;
  out << "\t\t\tEPS_SCF " << _settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) << std::endl;
  out << "\t\t\tMAX_SCF " << _settings.getInt(Utils::SettingsNames::maxScfIterations) << std::endl;
  out << "\t\t\tADDED_MOS " << _settings.getInt(SettingsNames::additionalMos) << std::endl;
  // subsections
  auto mixing = _settings.getString(Utils::SettingsNames::scfDamping);
  std::for_each(mixing.begin(), mixing.end(), [](char& c) { c = ::toupper(c); });
  if (!mixing.empty() && mixing != "NONE") {
    out << "\t\t\t&MIXING T" << std::endl;
    out << "\t\t\t\tMETHOD " << mixing << std::endl;
    out << "\t\t\t&END MIXING" << std::endl;
  }
  auto eTemp = _settings.getDouble(Utils::SettingsNames::electronicTemperature);
  if (eTemp > 0.0) {
    out << "\t\t\t&SMEAR ON" << std::endl;
    out << "\t\t\t\tMETHOD FERMI_DIRAC" << std::endl;
    out << "\t\t\t\tELECTRONIC_TEMPERATURE [K] " << eTemp << std::endl;
    out << "\t\t\t&END SMEAR" << std::endl;
  }
  auto ot = _settings.getString(SettingsNames::orbitalTransformation);
  if (!ot.empty()) {
    out << "\t\t\t&OT" << std::endl;
    out << "\t\t\t\tMINIMIZER " << ot << std::endl;
    out << "\t\t\t\tPRECONDITIONER FULL_ALL" << std::endl;
    out << "\t\t\t&END OT" << std::endl;
  }
  int outerScf = _settings.getInt(SettingsNames::outerScf);
  if (outerScf > 0) {
    out << "\t\t\t&OUTER_SCF" << std::endl;
    out << "\t\t\t\tMAX_SCF " << outerScf << std::endl;
    out << "\t\t\t\tEPS_SCF " << _settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) << std::endl;
    out << "\t\t\t&END OUTER_SCF" << std::endl;
  }
  out << "\t\t&END SCF" << std::endl;
}

void Cp2kInputFileCreator::printPoissonSolver(std::ostream& out) const {
  auto psolver = _settings.getString(SettingsNames::poissonSolver);
  if (!psolver.empty()) {
    out << "\t\t&POISSON" << std::endl;
    out << "\t\t\tPSOLVER " << psolver << std::endl;
    out << "\t\t&END POISSON" << std::endl;
  }
}

void Cp2kInputFileCreator::printGridInput(std::ostream& out) const {
  out << "\t\t&MGRID" << std::endl;
  out << "\t\t\tNGRIDS " << _settings.getInt(SettingsNames::nGrids) << std::endl;
  out << "\t\t\tCUTOFF " << _settings.getDouble(SettingsNames::planeWaveCutoff) << std::endl;
  out << "\t\t\tREL_CUTOFF " << _settings.getDouble(SettingsNames::relMultiGridCutoff) << std::endl;
  out << "\t\t&END MGRID" << std::endl;
}

void Cp2kInputFileCreator::printMatrixPrint(std::ostream& out) const {
  if (_requiredProperties.containsSubSet(Property::DensityMatrix) ||
      _requiredProperties.containsSubSet(Property::OverlapMatrix) ||
      _requiredProperties.containsSubSet(Property::BondOrderMatrix)) {
    out << "\t\t&PRINT" << std::endl;
    out << "\t\t\t&AO_MATRICES" << std::endl;
    auto addOutput = _settings.getString(SettingsNames::additionalOutputFile);
    if (!addOutput.empty() && addOutput != _settings.getString(SettingsNames::cp2kFilenameBase)) {
      out << "\t\t\t\tFILENAME ./" << addOutput << std::endl;
    }
    out << "\t\t\t\tDENSITY" << std::endl;
    out << "\t\t\t\tOVERLAP" << std::endl;
    out << "\t\t\t&END AO_MATRICES" << std::endl;
    out << "\t\t&END PRINT" << std::endl;
  }
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
