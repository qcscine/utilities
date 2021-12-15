/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Cp2kInputFileCreator.h"
#include "Cp2kCalculatorSettings.h"
#include <Utils/CalculatorBasics/PropertyList.h>
#include <Utils/DataStructures/PeriodicBoundaries.h>
#include <Utils/Geometry/AtomCollection.h>
#include <Utils/IO/MolecularTrajectoryIO.h>
#include <Utils/Scf/LcaoUtils/SpinMode.h>
#include <algorithm>
#include <fstream>
#include <iterator>

namespace Scine {
namespace Utils {
namespace ExternalQC {

void Cp2kInputFileCreator::createInputFile(const std::string& filename, const AtomCollection& atoms, const Settings& settings,
                                           const PropertyList& requiredProperties, const std::string& basename) const {
  std::ofstream fout;
  fout.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  fout.open(filename);
  createInputFile(fout, atoms, settings, requiredProperties, basename);
  fout.close();
}

void Cp2kInputFileCreator::createInputFile(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                           const PropertyList& requiredProperties, const std::string& basename) const {
  printGlobal(out, requiredProperties, basename);
  printForceEval(out, atoms, settings, requiredProperties);
}

void Cp2kInputFileCreator::printGlobal(std::ostream& out, const PropertyList& requiredProperties,
                                       const std::string& basename) const {
  std::string runType = (requiredProperties.containsSubSet(Property::Hessian)) ? "VIBRATIONAL_ANALYSIS" : "ENERGY_FORCE";
  std::string printLevel = "MEDIUM";
  out << "&GLOBAL" << std::endl;
  out << "\tPROJECT " << basename << std::endl;
  out << "\tRUN_TYPE " << runType << std::endl;
  out << "\tPRINT_LEVEL " << printLevel << std::endl;
  out << "\tEXTENDED_FFT_LENGTHS" << std::endl; // necessary for slab calculations
  out << "\t&PRINT" << std::endl;
  out << "\t\tPHYSCON FALSE" << std::endl; // no physical constants
  out << "\t&END PRINT" << std::endl;
  out << "&END GLOBAL" << std::endl << std::endl;
}

void Cp2kInputFileCreator::printForceEval(std::ostream& out, const AtomCollection& atoms, const Settings& settings,
                                          const PropertyList& requiredProperties) const {
  out << "&FORCE_EVAL" << std::endl;
  if (requiredProperties.containsSubSet(Property::Gradients)) {
    out << "\t&PRINT" << std::endl;
    out << "\t\t&FORCES ON" << std::endl;
    out << "\t\t\tNDIGITS 18" << std::endl;
    out << "\t\t&END FORCES" << std::endl;
    out << "\t&END PRINT" << std::endl;
  }
  out << "\tMETHOD QUICKSTEP" << std::endl;
  printSubsys(out, atoms, settings);
  // add method variations here if something other than DFT has been implemented
  printDftInput(out, settings, requiredProperties);
  out << "&END FORCE_EVAL" << std::endl;
}

void Cp2kInputFileCreator::printSubsys(std::ostream& out, const AtomCollection& atoms, const Settings& settings) const {
  out << "\t&SUBSYS" << std::endl;
  out << "\t\t&PRINT" << std::endl;
  out << "\t\t\t&SYMMETRY" << std::endl;
  out << "\t\t\t\tALL" << std::endl;
  if (settings.getString(Utils::SettingsNames::periodicBoundaries) == SettingsNames::moleculePeriodicBoundaries) {
    out << "\t\t\t\tMOLECULE" << std::endl;
  }
  out << "\t\t\t&END SYMMETRY" << std::endl;
  out << "\t\t&END PRINT" << std::endl;
  printCell(out, settings);
  printCoords(out, atoms);
  printBasis(out, atoms, settings);
  out << "\t&END SUBSYS" << std::endl;
}

void Cp2kInputFileCreator::printCell(std::ostream& out, const Settings& settings) const {
  out << "\t\t&CELL" << std::endl;
  auto pbc = PeriodicBoundaries(settings.getString(Utils::SettingsNames::periodicBoundaries));
  out << "\t\t\tABC " << std::to_string(pbc.getA().norm() * Constants::angstrom_per_bohr) << " "
      << std::to_string(pbc.getB().norm() * Constants::angstrom_per_bohr) << " "
      << std::to_string(pbc.getC().norm() * Constants::angstrom_per_bohr) << std::endl;
  out << "\t\t\tALPHA_BETA_GAMMA " << std::to_string(pbc.getAlpha()) << " " << std::to_string(pbc.getBeta()) << " "
      << std::to_string(pbc.getGamma()) << std::endl;
  out << "\t\t\tPERIODIC XYZ" << std::endl;
  out << "\t\t&END CELL" << std::endl;
}

void Cp2kInputFileCreator::printCoords(std::ostream& out, const AtomCollection& atoms) const {
  out << "\t\t&COORD" << std::endl;
  for (const auto& a : atoms) {
    out << "\t\t";
    MolecularTrajectoryIO::writeXYZLine(out, a.getElementType(), a.getPosition());
  }
  out << "\t\t&END COORD" << std::endl;
  out << "\t\t&TOPOLOGY\n\t\t\t&CENTER_COORDINATES\n\t\t\t&END\n\t\t&END TOPOLOGY" << std::endl; // centering coords
}

void Cp2kInputFileCreator::printBasis(std::ostream& out, const AtomCollection& atoms, const Settings& settings) const {
  auto basisEntry = settings.getString(Utils::SettingsNames::basisSet);
  std::for_each(basisEntry.begin(), basisEntry.end(), [](char& c) { c = ::toupper(c); });
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
  // get unique elements
  auto uniqueElements = atoms.getElements();
  std::sort(uniqueElements.begin(), uniqueElements.end());
  ElementTypeCollection::iterator it;
  it = std::unique(uniqueElements.begin(), uniqueElements.end());
  uniqueElements.resize(std::distance(uniqueElements.begin(), it));
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
  // write different elements
  for (const auto& element : uniqueElements) {
    std::string actualBasis = zetaInfo + "-MOLOPT";
    if (std::find(_highQualityBasisElements.begin(), _highQualityBasisElements.end(), element) ==
        _highQualityBasisElements.end()) {
      actualBasis += "-SR";
    }
    actualBasis += "-GTH";
    auto method = settings.getString(Utils::SettingsNames::method);
    std::for_each(method.begin(), method.end(), [](char& c) { c = ::toupper(c); });

    out << "\t\t&KIND " << ElementInfo::symbol(element) << std::endl;
    out << "\t\t\tELEMENT " << ElementInfo::symbol(element) << std::endl;
    out << "\t\t\tBASIS_SET " << actualBasis << std::endl;
    out << "\t\t\tPOTENTIAL GTH-" + method << std::endl;
    out << "\t\t&END KIND" << std::endl;
  }
}

void Cp2kInputFileCreator::printDftInput(std::ostream& out, const Settings& settings,
                                         const PropertyList& requiredProperties) const {
  out << "\t&DFT" << std::endl;
  printElectronicStructureBasics(out, settings);
  printScfInput(out, settings);
  printPoissonSolver(out, settings);
  printGridInput(out, settings);
  printMatrixPrint(out, settings, requiredProperties);
  out << "\t&END DFT" << std::endl;
}

void Cp2kInputFileCreator::printElectronicStructureBasics(std::ostream& out, const Settings& settings) const {
  out << "\t\tCHARGE " << settings.getInt(Utils::SettingsNames::molecularCharge) << std::endl;
  int multiplicity = settings.getInt(Utils::SettingsNames::spinMultiplicity);
  out << "\t\tMULTIPLICITY " << multiplicity << std::endl;
  auto spinMode = SpinModeInterpreter::getSpinModeFromString(settings.getString(Utils::SettingsNames::spinMode));
  out << "\t\t" << determineCp2kSpinMode(spinMode, multiplicity) << std::endl;
  out << "\t\tBASIS_SET_FILE_NAME BASIS_MOLOPT" << std::endl; // TODO modify if more basis sets are implemented
  out << "\t\t&XC" << std::endl;
  out << "\t\t\t&XC_FUNCTIONAL " << settings.getString(Utils::SettingsNames::method) << std::endl;
  out << "\t\t\t&END XC_FUNCTIONAL" << std::endl;
  printDispersionCorrection(out, settings);
  if (settings.getString(Utils::SettingsNames::periodicBoundaries) == SettingsNames::moleculePeriodicBoundaries) {
    out << "\t\t\t&XC_GRID" << std::endl;
    out << "\t\t\t\tXC_SMOOTH_RHO NN50" << std::endl;
    out << "\t\t\t\tXC_DERIV NN50_SMOOTH" << std::endl;
    out << "\t\t\t&END XC_GRID" << std::endl;
  }
  out << "\t\t&END XC" << std::endl;
  if (settings.getBool(SettingsNames::dipoleCorrection)) {
    out << "\t\tSURFACE_DIPOLE_CORRECTION" << std::endl;
  }
}

void Cp2kInputFileCreator::printDispersionCorrection(std::ostream& out, const Settings& settings) const {
  auto vdwFunc = settings.getString(SettingsNames::vdwFunctional);
  if (!vdwFunc.empty()) {
    out << "\t\t\t&VDW_POTENTIAL" << std::endl;
    if (vdwFunc.find("DFTD") != std::string::npos) {
      out << "\t\t\t\tPOTENTIAL_TYPE PAIR_POTENTIAL" << std::endl;
      out << "\t\t\t\t&PAIR_POTENTIAL" << std::endl;
      out << "\t\t\t\t\tTYPE " << vdwFunc << std::endl;
      out << "\t\t\t\t\tREFERENCE_FUNCTIONAL " << settings.getString(Utils::SettingsNames::method) << std::endl;
      out << "\t\t\t\t\tPARAMETER_FILE_NAME dftd3.dat" << std::endl;
      out << "\t\t\t\t&END PAIR_POTENTIAL" << std::endl;
    }
    else {
      out << "\t\t\t\tPOTENTIAL_TYPE NON_LOCAL" << std::endl;
      out << "\t\t\t\t&NON_LOCAL" << std::endl;
      out << "\t\t\t\t\tTYPE " << vdwFunc << std::endl;
      out << "\t\t\t\t\tKERNEL_FILE_NAME vdW_kernel_table.dat" << std::endl;
      out << "\t\t\t\t&END NON_LOCAL" << std::endl;
    }
    out << "\t\t\t&END VDW_POTENTIAL" << std::endl;
  }
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
                         " in settings."); // this should have been handled by settings
}

void Cp2kInputFileCreator::printScfInput(std::ostream& out, const Settings& settings) const {
  // E accuracy
  out << "\t\t&SCF" << std::endl;
  // single keywords
  out << "\t\t\tSCF_GUESS " << settings.getString(SettingsNames::scfGuess) << std::endl;
  out << "\t\t\tEPS_SCF " << settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) << std::endl;
  out << "\t\t\tMAX_SCF " << settings.getInt(Utils::SettingsNames::maxScfIterations) << std::endl;
  out << "\t\t\tADDED_MOS " << settings.getInt(SettingsNames::additionalMos) << std::endl;
  // subsections
  auto mixing = settings.getString(Utils::SettingsNames::scfDamping);
  std::for_each(mixing.begin(), mixing.end(), [](char& c) { c = ::toupper(c); });
  if (!mixing.empty() && mixing != "NONE") {
    out << "\t\t\t&MIXING T" << std::endl;
    out << "\t\t\t\tMETHOD " << mixing << std::endl;
    out << "\t\t\t&END MIXING" << std::endl;
  }
  auto eTemp = settings.getDouble(Utils::SettingsNames::electronicTemperature);
  if (eTemp > 0.0) {
    out << "\t\t\t&SMEAR ON" << std::endl;
    out << "\t\t\t\tMETHOD FERMI_DIRAC" << std::endl;
    out << "\t\t\t\tELECTRONIC_TEMPERATURE [K] " << eTemp << std::endl;
    out << "\t\t\t&END SMEAR" << std::endl;
  }
  auto ot = settings.getString(SettingsNames::orbitalTransformation);
  if (!ot.empty()) {
    out << "\t\t\t&OT" << std::endl;
    out << "\t\t\t\tMINIMIZER " << ot << std::endl;
    out << "\t\t\t\tPRECONDITIONER FULL_ALL" << std::endl;
    out << "\t\t\t&END OT" << std::endl;
  }
  int outerScf = settings.getInt(SettingsNames::outerScf);
  if (outerScf > 0) {
    out << "\t\t\t&OUTER_SCF" << std::endl;
    out << "\t\t\t\tMAX_SCF " << outerScf << std::endl;
    out << "\t\t\t\tEPS_SCF " << settings.getDouble(Utils::SettingsNames::selfConsistenceCriterion) << std::endl;
    out << "\t\t\t&END OUTER_SCF" << std::endl;
  }
  out << "\t\t&END SCF" << std::endl;
}

void Cp2kInputFileCreator::printPoissonSolver(std::ostream& out, const Settings& settings) const {
  auto psolver = settings.getString(SettingsNames::poissonSolver);
  if (!psolver.empty()) {
    out << "\t\t&POISSON" << std::endl;
    out << "\t\t\tPSOLVER " << psolver << std::endl;
    out << "\t\t&END POISSON" << std::endl;
  }
}

void Cp2kInputFileCreator::printGridInput(std::ostream& out, const Settings& settings) const {
  out << "\t\t&MGRID" << std::endl;
  out << "\t\t\tNGRIDS " << settings.getInt(SettingsNames::nGrids) << std::endl;
  out << "\t\t\tCUTOFF " << settings.getDouble(SettingsNames::planeWaveCutoff) << std::endl;
  out << "\t\t\tREL_CUTOFF " << settings.getDouble(SettingsNames::relMultiGridCutoff) << std::endl;
  out << "\t\t&END MGRID" << std::endl;
}

void Cp2kInputFileCreator::printMatrixPrint(std::ostream& out, const Settings& settings,
                                            const PropertyList& requiredProperties) const {
  if (requiredProperties.containsSubSet(Property::DensityMatrix) || requiredProperties.containsSubSet(Property::OverlapMatrix) ||
      requiredProperties.containsSubSet(Property::BondOrderMatrix)) {
    out << "\t\t&PRINT" << std::endl;
    out << "\t\t\t&AO_MATRICES" << std::endl;
    auto addOutput = settings.getString(SettingsNames::additionalOutputFile);
    if (!addOutput.empty() && addOutput != settings.getString(SettingsNames::cp2kFilenameBase)) {
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
