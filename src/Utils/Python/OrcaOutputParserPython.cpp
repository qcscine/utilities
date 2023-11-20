/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include <Utils/Bonds/BondOrderCollection.h>
#include <Utils/DataStructures/SingleParticleEnergies.h>
#include <Utils/ExternalQC/Orca/OrcaMainOutputParser.h>
#include <Utils/Pybind.h>

void init_orca_parser(pybind11::module& m) {
  using namespace Scine::Utils;
  using namespace ExternalQC;
  using Parser = ExternalQC::OrcaMainOutputParser;

  pybind11::class_<Parser> parser(m, "OrcaOutputParser");
  parser.def(pybind11::init<std::string>(), pybind11::arg("filename"));
  parser.def("raise_errors", &Parser::checkForErrors, "Raise an exception for any errors");
  parser.def("energy", &Parser::getEnergy, "Parse energy in Hartree from output");
  parser.def("bond_orders", &Parser::getBondOrders, "Parse Mayber bond orders from output");
  parser.def("hirshfeld_charges", &Parser::getHirshfeldCharges, "Parse Hirshfeld charges from output");
  parser.def("num_atoms", &Parser::getNumberAtoms, "Parse number of atoms from output");
  parser.def("gradients", &Parser::getGradients, "Parse gradients from output");
  parser.def("orbital_energies", &Parser::getOrbitalEnergies, "Parse orbital energies in Hartree from output");
  parser.def("temperature", &Parser::getTemperature, "Parse temperature in Kelvin from output");
  parser.def("enthalpy", &Parser::getEnthalpy, "Parse enthalpy in Hartree from output");
  parser.def("entropy", &Parser::getEntropy, "Parse entropy in Hartree/Kelvin from output");
  parser.def("zpve", &Parser::getZeroPointVibrationalEnergy, "Parse zero point vibrational energy in Hartree from output");
  parser.def("symmetry_number", &Parser::getSymmetryNumber,
             "Parse the molecular symmetry number for which thermochemistry was computed from the output");
  parser.def("moessbauer_asymmetry_parameter", &Parser::getMoessbauerAsymmetryParameter, "");
  parser.def("moessbauer_quadrupole_splittings", &Parser::getMoessbauerQuadrupoleSplittings, "");
  parser.def("moessbauer_iron_electron_densities", &Parser::getMoessbauerIronElectronDensities, "");
}
