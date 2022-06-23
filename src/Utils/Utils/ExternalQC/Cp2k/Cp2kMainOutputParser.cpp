/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Cp2kMainOutputParser.h"
#include <Utils/Constants.h>
#include <Utils/ExternalQC/Exceptions.h>
#include <Utils/Geometry/ElementInfo.h>
#include <Utils/IO/Regex.h>
#include <Utils/Scf/LcaoUtils/LcaoUtils.h>
#include <fstream>
#include <numeric>

namespace Scine {
namespace Utils {
namespace ExternalQC {

Cp2kMainOutputParser::Cp2kMainOutputParser(const std::string& filename, const std::string& additionalOutput) {
  content_ = extractContent(filename);
  if (!additionalOutput.empty()) {
    additionalContent_ = extractContent(additionalOutput);
  }
  extractRuntype();
}

std::string Cp2kMainOutputParser::extractContent(const std::string& filename) {
  std::ifstream fin;
  fin.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  // TODO: This file is never properly closed
  fin.open(filename);
  return std::string(std::istreambuf_iterator<char>{fin}, {});
}

void Cp2kMainOutputParser::extractRuntype() {
  std::regex rgx("Run type\\s+(\\w+)\\s");
  std::smatch match;
  std::sregex_iterator iter(content_.begin(), content_.end(), rgx);
  std::sregex_iterator end;

  if (iter == end || iter->size() != 2) {
    throw OutputFileParsingError("Could not read run type of calculation.");
  }
  runtype_ = (*iter)[1];
}

void Cp2kMainOutputParser::checkForErrors() const {
  std::regex regex("WARNING.+SCF run NOT converged");
  std::smatch match;
  if (std::regex_search(content_, match, regex)) {
    throw ScfNotConvergedError("SCF in CP2K calculation did not converge.");
  }
  std::regex regex2(R"(E\s?R\s?R\s?O\s?R\s?)");
  std::smatch match2;
  if (std::regex_search(content_, match2, regex2)) {
    throw OutputFileParsingError("CP2K encountered an error during the calculation.");
  }
}

double Cp2kMainOutputParser::getEnergy() const {
  std::string energyRegexString = "ENERGY. Total FORCE_EVAL \\( QS \\) energy .a\\.u\\..\\:\\s+";
  std::string hessianRegexString = "Minimum Structure - Energy and Forces:\\s+VIB.\\s+Total Energy:\\s+";
  std::string searchString = (runtype_ == "VIBRATIONAL_ANALYSIS") ? hessianRegexString : energyRegexString;
  std::regex regex(searchString + Regex::capturingFloatingPointNumber());
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b) {
    throw OutputFileParsingError("Energy could not be read from CP2K output.");
  }

  return std::stod(matches[1]);
}

GradientCollection Cp2kMainOutputParser::getGradients() const {
  std::string energyRegexString = "ATOMIC FORCES in .a\\.u\\..\\s+#\\s+Atom\\s+Kind\\s+Element\\s+X\\s+Y\\s+Z\\s+"
                                  "((?:\\d+\\s+\\d+\\s+" +
                                  Regex::elementSymbol() + "\\s+" + Regex::floatingPointNumber() + "\\s+" +
                                  Regex::floatingPointNumber() + "\\s+" + Regex::floatingPointNumber() + "\\s+" +
                                  ")+)SUM OF ATOMIC FORCES";
  std::string hessianRegexString =
      "Minimum Structure - Energy and Forces:\\s+VIB.\\s+Total Energy:\\s+" + Regex::floatingPointNumber() +
      "\\s+"
      "VIB.\\s+ATOM\\s+X\\s+Y\\s+Z\\s+"
      "((?:VIB.\\s+" +
      Regex::elementSymbol() + "\\s+" + Regex::floatingPointNumber() + "\\s+" + Regex::floatingPointNumber() + "\\s+" +
      Regex::floatingPointNumber() + "\\s+" + ")+)VIB.\\s+Hessian in cartesian coordinates";
  std::string searchString = (runtype_ == "VIBRATIONAL_ANALYSIS") ? hessianRegexString : energyRegexString;
  // capture correct forces block
  std::regex regex(searchString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b) {
    throw OutputFileParsingError("Gradients could not be read from CP2K output.");
  }
  std::string forcesBlock = matches[1];
  // capture all the numbers
  searchString = Regex::capturingFloatingPointNumber() + "\\s+" + Regex::capturingFloatingPointNumber() + "\\s+" +
                 Regex::capturingFloatingPointNumber() + "\\s+";
  std::regex regex2(searchString);
  std::sregex_iterator iter(forcesBlock.begin(), forcesBlock.end(), regex2);
  std::sregex_iterator end;
  std::vector<double> numbers;
  while (iter != end) {
    if (iter->size() != 4) {
      throw OutputFileParsingError("Gradients could not be read from CP2K output.");
    }
    for (unsigned i = 1; i < 4; ++i) {
      numbers.push_back(std::stod((*iter)[i]));
    }
    ++iter;
  }
  GradientCollection result = Eigen::Map<GradientCollection>(numbers.data(), numbers.size() / 3, 3);
  result *= -1.0;

  return result;
}

HessianMatrix Cp2kMainOutputParser::getHessian() const {
  // get n atoms
  std::regex regex("Atomic kind:\\s+" + Regex::elementSymbol() + "\\s+Number of atoms:\\s+" + Regex::capturingIntegerNumber());
  std::sregex_iterator iter(content_.begin(), content_.end(), regex);
  std::sregex_iterator end;
  std::vector<int> numbers;
  while (iter != end) {
    if (iter->size() != 2) {
      throw OutputFileParsingError("Gradients could not be read from CP2K output.");
    }
    numbers.push_back(std::stoi((*iter)[1]));
    ++iter;
  }
  int nAtoms = std::accumulate(numbers.begin(), numbers.end(), decltype(numbers)::value_type(0));
  // find Hessian block
  std::string endKeyword = "VIB";
  // strings declarations necessary for in-line comments
  std::string searchString = std::string("Hessian in cartesian coordinates") + // beginning identifier
                             std::string("((?:") +                             // begin multiple blocks
                             std::string("(?:\\s+\\d+)+\\s+") + // column declarations + row number in next line
                             // element symbol and one or more numbers
                             std::string(Regex::elementSymbol() + "(?:\\s+" + Regex::floatingPointNumber() + ")+") +
                             std::string(")+)") +              // end multiple blocks
                             std::string("\\s+" + endKeyword); // end identifier
  std::regex regex2(searchString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex2);
  if (!b) {
    throw OutputFileParsingError("Hessian could not be read from CP2K output.");
  }
  std::string hessianBlock = matches[1];

  auto hessian = parseMatrixFromStringBlock(hessianBlock, "Hessian", nAtoms * 3);
  if (hessian.isApprox(HessianMatrix::Zero(nAtoms * 3, nAtoms * 3))) {
    throw OutputFileParsingError("Hessian could not be read from CP2K output.");
  }
  return hessian;
}

std::vector<double> Cp2kMainOutputParser::getHirshfeldCharges() const {
  // parse output line by line
  std::regex beginRegex("Hirshfeld Charges");
  std::regex endRegex("Total Charge");
  std::regex numberRegex(Regex::floatingPointNumber());
  std::stringstream stream(content_);
  std::string line;
  std::vector<double> charges;
  bool parsing = false;
  while (std::getline(stream, line)) {
    if (!line.empty()) {
      std::smatch match;
      if (std::regex_search(line, match, beginRegex)) { // we have to start parsing
        if (parsing || !charges.empty()) {              // we already parsed and got beginning again, something is weird
          throw OutputFileParsingError("Charges could not be read from CP2K output.");
        }
        parsing = true;
      }
      else if (std::regex_search(line, match, endRegex)) { // we have to end parsing
        if (!parsing || charges.empty()) {                 // we did not parse and got end, something is weird
          throw OutputFileParsingError("Charges could not be read from CP2K output.");
        }
        parsing = false;
      }
      else if (parsing && std::regex_search(line, match, numberRegex)) { // we have a line specifying a charge
        double currentCharge = 0.0;
        // capture all floats in line and pick last one
        std::regex numberRegexCapture(Regex::capturingFloatingPointNumber());
        std::sregex_iterator it(line.begin(), line.end(), numberRegexCapture);
        std::sregex_iterator end;
        while (it != end) {
          if (it->size() != 2) {
            throw OutputFileParsingError("Charges could not be read from CP2K output.");
          }
          currentCharge = std::stod((*it)[1]);
          ++it;
        }
        charges.push_back(currentCharge); // last number in line is partial charge
      }
    }
  }
  if (charges.empty()) { // we did not parse and got end, something is weird
    throw OutputFileParsingError("Charges could not be read from CP2K output.");
  }
  return charges;
}

std::vector<int> Cp2kMainOutputParser::getGridCounts() const {
  std::string searchstring = "count for grid\\s+\\d+:\\s+" + Regex::capturingIntegerNumber() +
                             "\\s+cutoff .a\\.u\\..\\s+" + Regex::floatingPointNumber();
  std::regex regex(searchstring);
  std::sregex_iterator iter(content_.begin(), content_.end(), regex);
  std::sregex_iterator end;
  std::vector<int> gridCounts;
  while (iter != end) {
    if (iter->size() != 2) {
      throw OutputFileParsingError("Grid counts could not be read from CP2K output.");
    }
    gridCounts.push_back(std::stoi((*iter)[1]));
    ++iter;
  }
  return gridCounts;
}

BondOrderCollection Cp2kMainOutputParser::getBondOrders(const ElementTypeCollection& elements, SpinMode spinMode) const {
  try {
    auto density = getDensityMatrix(spinMode);
    auto indices = getAtomAoIndex(elements);
    auto overlap = getOverlapMatrix();
    auto boCollection = BondOrderCollection(elements.size());
    LcaoUtils::calculateBondOrderMatrix(boCollection, density, overlap, indices);
    return boCollection;
  }
  catch (const OutputFileParsingError& error) {
    std::string underlyingError = error.what();
    throw OutputFileParsingError("Bond orders could not be read from CP2K output, because " + underlyingError);
  }
}

DensityMatrix Cp2kMainOutputParser::getDensityMatrix(SpinMode spinMode) const {
  int nAOs = getNumberOfAos();
  auto nElectrons = getNumberOfElectrons();
  // either one block (restricted) or 2 blocks for alpha and beta
  // therefore this functions always uses vectors with either just 1 entry or 2 entries
  std::vector<std::string> restrictedHeaders = {"DENSITY MATRIX"};
  std::vector<std::string> unrestrictedHeaders = {"DENSITY MATRIX FOR ALPHA SPIN", "DENSITY MATRIX FOR BETA SPIN"};
  std::vector<std::string> searchHeaders = (spinMode == SpinMode::Restricted) ? restrictedHeaders : unrestrictedHeaders;
  std::vector<std::string> densityBlocks;
  std::unique_ptr<std::string> content = (additionalContent_.empty()) ? std::make_unique<std::string>(content_)
                                                                      : std::make_unique<std::string>(additionalContent_);
  for (const auto& header : searchHeaders) {
    auto densityBlock = parseMatrixStringBlockFromFileContent(*content, std::regex(header));
    if (densityBlock.empty()) {
      throw OutputFileParsingError("Density matrix could not be read from CP2K output.");
    }
    densityBlocks.push_back(densityBlock);
  }
  // extract matrices from block
  std::vector<Eigen::MatrixXd> matrices;
  for (const auto& densityBlock : densityBlocks) {
    auto densityMatrix = parseMatrixFromStringBlock(densityBlock, "Density Matrix", nAOs);
    matrices.push_back(densityMatrix);
  }
  // translate Eigen Matrix to DensityMatrix datastructure
  auto densityMatrix = DensityMatrix();
  if (spinMode == SpinMode::Restricted) {
    if (matrices.size() != 1) {
      throw OutputFileParsingError("Density matrix could not be read from CP2K output.");
    }
    densityMatrix.setDensity(std::move(matrices[0]), std::move(nElectrons[0]));
  }
  else {
    if (matrices.size() != 2) {
      throw OutputFileParsingError("Density matrix could not be read from CP2K output.");
    }
    densityMatrix.setDensity(std::move(matrices[0]), std::move(matrices[1]), std::move(nElectrons[0]),
                             std::move(nElectrons[1]));
  }
  return densityMatrix;
}

Eigen::MatrixXd Cp2kMainOutputParser::getOverlapMatrix() const {
  int nAOs = getNumberOfAos();
  std::string searchString = "OVERLAP MATRIX";
  std::unique_ptr<std::string> content = (additionalContent_.empty()) ? std::make_unique<std::string>(content_)
                                                                      : std::make_unique<std::string>(additionalContent_);
  auto overlapBlock = parseMatrixStringBlockFromFileContent(*content, std::regex(searchString));
  if (overlapBlock.empty()) {
    throw OutputFileParsingError("Overlap matrix could not be read from CP2K output.");
  }
  return parseMatrixFromStringBlock(overlapBlock, "Overlap matrix", nAOs);
}

AtomsOrbitalsIndexes Cp2kMainOutputParser::getAtomAoIndex(const ElementTypeCollection& elements) const {
  std::string blockSearchString = "Atomic kind:\\s+" + Regex::capturingElementSymbol() +
                                  "\\s+Number of atoms:\\s+"
                                  "\\d+\\s+.+\\s+(?:(?:\\s+\\w+)+:\\s+\\d+)+";
  std::string sphericalsSearchString = "Number of spherical basis functions:\\s+" + Regex::capturingIntegerNumber();
  std::regex regex(blockSearchString);
  std::regex sphericalsRegex(sphericalsSearchString);
  std::sregex_iterator iter(content_.begin(), content_.end(), regex);
  std::sregex_iterator end;
  std::map<Utils::ElementType, int> nSphericals;
  // cycle through element blocks
  while (iter != end) {
    if (iter->size() != 2) {
      throw OutputFileParsingError("AO indices could not be read from CP2K output.");
    }
    Utils::ElementType elementType = Utils::ElementInfo::elementTypeForSymbol((*iter)[1]);
    // search block to get the number of sphericals
    std::string block = (*iter)[0];
    std::smatch matches;
    bool b = std::regex_search(block, matches, sphericalsRegex);
    if (!b || matches.size() != 2) {
      throw OutputFileParsingError("AO indices could not be read from CP2K output.");
    }
    nSphericals.insert(std::make_pair(elementType, std::stoi(matches[1])));
    ++iter;
  }
  // check if all elements were found in output
  for (const auto& element : elements) {
    if (nSphericals.find(element) == nSphericals.end()) {
      throw OutputFileParsingError("AO indices could not be read from CP2K output.");
    }
  }
  auto result = AtomsOrbitalsIndexes(elements.size());
  for (const auto& element : elements) {
    result.addAtom(nSphericals.at(element));
  }
  assert(result.getNAtomicOrbitals() == getNumberOfAos());
  return result;
}

int Cp2kMainOutputParser::getNumberOfAos() const {
  std::string searchString = "Spherical basis functions:\\s+" + Regex::capturingIntegerNumber();
  std::regex regex(searchString);
  std::smatch matches;
  bool b = std::regex_search(content_, matches, regex);
  if (!b || matches.size() != 2) {
    throw OutputFileParsingError("Number of AOs could not be read from CP2K output.");
  }
  return std::stoi(matches[1]);
}

std::vector<int> Cp2kMainOutputParser::getNumberOfElectrons() const {
  // find out number of electrons
  std::string searchString = "Number of electrons:\\s+" + Regex::capturingIntegerNumber();
  std::regex regex(searchString);
  std::sregex_iterator iter(content_.begin(), content_.end(), regex);
  std::sregex_iterator end;
  std::vector<int> nElectrons;
  while (iter != end) {
    if (iter->size() != 2) {
      throw OutputFileParsingError("Number of electrons could not be read from CP2K output.");
    }
    nElectrons.push_back(std::stoi((*iter)[1]));
    ++iter;
  }
  return nElectrons;
}

Eigen::Matrix3d Cp2kMainOutputParser::getStressTensor() const {
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
  int rowCount = 0;
  std::regex startRegex("STRESS TENSOR");
  std::stringstream stream(content_);
  std::string line;
  bool parseStress = false;
  // parse output line by line
  while (std::getline(stream, line)) {
    if (!line.empty()) {
      std::smatch match;
      if (parseStress) { // we are reading matrix and getting the first 3 lines with 3 floating point numbers
        std::regex numbers("\\s+" + Regex::capturingFloatingPointNumber() + "\\s+" +
                           Regex::capturingFloatingPointNumber() + "\\s+" + Regex::capturingFloatingPointNumber());
        if (std::regex_search(line, match, numbers)) {
          if (match.size() != 4) {
            throw OutputFileParsingError("Stress tensor could not be read.");
          }
          for (int i = 0; i < 3; ++i) {
            stress(rowCount, i) = std::stod(match[i + 1]);
          }
          rowCount++;
          if (rowCount == 3) {
            // CP2K gives stress tensor in GPa
            return stress * 1e9 * Constants::hartree_per_joule * std::pow(Constants::meter_per_bohr, 3);
          }
        }
      }
      else if (std::regex_search(line, match, startRegex)) { // we have to start reading matrix
        parseStress = true;
      }
    }
  }
  throw OutputFileParsingError("Stress tensor could not be read.");
}

int Cp2kMainOutputParser::getSymmetryNumber() const {
  // TODO
  // CP2K does not return int number system and symmetry detection is not reliable
  return 1;
}

std::string Cp2kMainOutputParser::parseMatrixStringBlockFromFileContent(const std::string& content,
                                                                        const std::regex& startKeyword) const {
  // do not change this parsing to complete regex, unless you know what you are doing
  // std::regex_iterators might segfault for long files, which these matrix files can be
  // https://stackoverflow.com/questions/36304204/c-regex-segfault-on-long-sequences
  std::stringstream stream(content);
  std::string line;
  std::string currentMatrixBlock;
  int emptyLineCounter = 0;
  bool addingToBlock = false;
  // parse content line by line
  // general logic:
  // A matrix block is started with the given keyword, e.g., 'DENSITY MATRIX'
  // A matrix block ends if there are 2 consecutive empty lines
  // 1 empty line does not mean anything
  // We only want the very last matrix block in the file
  while (std::getline(stream, line)) {
    if (line.empty() && addingToBlock) {
      // might be end of block
      emptyLineCounter++;
      if (emptyLineCounter == 2) {
        // we have reached end
        addingToBlock = false;
      }
      continue;
    }
    if (!addingToBlock) {
      // check for start
      std::smatch match;
      if (std::regex_search(line, match, startKeyword)) {
        // we have a new block
        currentMatrixBlock = ""; // empty if already parsed a block before
        addingToBlock = true;
      }
      continue; // in any case we do not want this line in our block
    }
    // reached here if we are currently adding and do not have an empty line
    emptyLineCounter = 0; // reset in case we had an empty line before
    currentMatrixBlock += line;
    currentMatrixBlock += "\n";
  }
  return currentMatrixBlock;
}

Eigen::MatrixXd Cp2kMainOutputParser::parseMatrixFromStringBlock(const std::string& block, const std::string& propName,
                                                                 int dimension) const {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(dimension, dimension);
  std::regex elementRegex(Regex::elementSymbol());
  std::sregex_iterator end;
  std::stringstream stream(block);
  std::string line;
  std::vector<int> currentCols;
  // parse block line by line
  while (std::getline(stream, line)) {
    if (!line.empty()) {
      std::smatch match;
      if (!std::regex_search(line, match, elementRegex)) { // we have a line specifying columns
        currentCols.clear();
        std::regex number(Regex::capturingIntegerNumber());
        std::sregex_iterator it(line.begin(), line.end(), number);
        while (it != end) {
          if (it->size() != 2) {
            throw OutputFileParsingError(propName + " could not be read from CP2K output.");
          }
          currentCols.push_back(std::stoi((*it)[1]));
          ++it;
        }
      }
      else { // we have a line specifying matrix entries
        std::vector<double> currentEntries;
        std::regex number(Regex::capturingFloatingPointNumber());
        std::sregex_iterator it(line.begin(), line.end(), number);
        int row = 0;
        while (it != end) {
          if (it->size() != 2) {
            throw OutputFileParsingError(propName + " could not be read from CP2K output.");
          }
          if (row == 0) {
            row = std::stoi((*it)[1]); // captures integer as well, this is the current row
          }
          else {
            currentEntries.push_back(std::stod((*it)[1]));
          }
          ++it;
        }
        if (currentEntries.size() < currentCols.size()) { // sanity check
          throw OutputFileParsingError(propName + " could not be read from CP2K output.");
        }
        // write into matrix
        // entries can have false positive numbers in front of real entries -> take only the x last ones with x being
        // the number of current columns
        unsigned long shift = currentEntries.size() - currentCols.size();
        for (unsigned long i = 0; i < currentCols.size(); ++i) {
          // Cp2k output starts with 1
          matrix(row - 1, currentCols[i] - 1) = currentEntries[shift + i];
        }
      }
    }
  }
  return matrix;
}

} // namespace ExternalQC
} // namespace Utils
} // namespace Scine
