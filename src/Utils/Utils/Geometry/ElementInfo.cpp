/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/ElementData.h"

namespace Scine {
namespace Utils {

ElementType ElementInfo::elementTypeForSymbol(const std::string& symbol) {
  auto it = stringToElementType.find(symbol);
  if (it == stringToElementType.end())
    throw ElementSymbolNotFound(symbol);
  return it->second;
}

std::string ElementInfo::symbol(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].symbol();
}

double ElementInfo::mass(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].mass();
}

double ElementInfo::vdwRadius(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].vdWRadius();
}

int ElementInfo::Z(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].Z();
}

int ElementInfo::valElectrons(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].valElectrons();
}

int ElementInfo::sElectrons(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].sElectrons();
}

int ElementInfo::pElectrons(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].pElectrons();
}

int ElementInfo::dElectrons(ElementType e) {
  return Constants::ElementDataSingleton::instance()[e].dElectrons();
}

std::map<std::string, ElementType> ElementInfo::stringToElementType = {
    {"H", ElementType::H},       {"He", ElementType::He},    {"Li", ElementType::Li}, {"Be", ElementType::Be},
    {"B", ElementType::B},       {"C", ElementType::C},      {"N", ElementType::N},   {"O", ElementType::O},
    {"F", ElementType::F},       {"Ne", ElementType::Ne},    {"Na", ElementType::Na}, {"Mg", ElementType::Mg},
    {"Al", ElementType::Al},     {"Si", ElementType::Si},    {"P", ElementType::P},   {"S", ElementType::S},
    {"Cl", ElementType::Cl},     {"Ar", ElementType::Ar},    {"K", ElementType::K},   {"Ca", ElementType::Ca},
    {"Sc", ElementType::Sc},     {"Ti", ElementType::Ti},    {"V", ElementType::V},   {"Cr", ElementType::Cr},
    {"Mn", ElementType::Mn},     {"Fe", ElementType::Fe},    {"Co", ElementType::Co}, {"Ni", ElementType::Ni},
    {"Cu", ElementType::Cu},     {"Zn", ElementType::Zn},    {"Ga", ElementType::Ga}, {"Ge", ElementType::Ge},
    {"As", ElementType::As},     {"Se", ElementType::Se},    {"Br", ElementType::Br}, {"Kr", ElementType::Kr},
    {"Rb", ElementType::Rb},     {"Sr", ElementType::Sr},    {"Y", ElementType::Y},   {"Zr", ElementType::Zr},
    {"Nb", ElementType::Nb},     {"Mo", ElementType::Mo},    {"Tc", ElementType::Tc}, {"Ru", ElementType::Ru},
    {"Rh", ElementType::Rh},     {"Pd", ElementType::Pd},    {"Ag", ElementType::Ag}, {"Cd", ElementType::Cd},
    {"In", ElementType::In},     {"Sn", ElementType::Sn},    {"Sb", ElementType::Sb}, {"Te", ElementType::Te},
    {"I", ElementType::I},       {"Xe", ElementType::Xe},    {"Cs", ElementType::Cs}, {"Ba", ElementType::Ba},
    {"La", ElementType::La},     {"Ce", ElementType::Ce},    {"Pr", ElementType::Pr}, {"Nd", ElementType::Nd},
    {"Pm", ElementType::Pm},     {"Sm", ElementType::Sm},    {"Eu", ElementType::Eu}, {"Gd", ElementType::Gd},
    {"Tb", ElementType::Tb},     {"Dy", ElementType::Dy},    {"Ho", ElementType::Ho}, {"Er", ElementType::Er},
    {"Tm", ElementType::Tm},     {"Yb", ElementType::Yb},    {"Lu", ElementType::Lu}, {"Hf", ElementType::Hf},
    {"Ta", ElementType::Ta},     {"W", ElementType::W},      {"Re", ElementType::Re}, {"Os", ElementType::Os},
    {"Ir", ElementType::Ir},     {"Pt", ElementType::Pt},    {"Au", ElementType::Au}, {"Hg", ElementType::Hg},
    {"Tl", ElementType::Tl},     {"Pb", ElementType::Pb},    {"Bi", ElementType::Bi}, {"Po", ElementType::Po},
    {"At", ElementType::At},     {"Rn", ElementType::Rn},    {"Fr", ElementType::Fr}, {"Ra", ElementType::Ra},
    {"Ac", ElementType::Ac},     {"Th", ElementType::Th},    {"Pa", ElementType::Pa}, {"U", ElementType::U},
    {"Np", ElementType::Np},     {"Pu", ElementType::Pu},    {"Am", ElementType::Am}, {"Cm", ElementType::Cm},
    {"Bk", ElementType::Bk},     {"Cf", ElementType::Cf},    {"Es", ElementType::Es}, {"Fm", ElementType::Fm},
    {"Md", ElementType::Md},     {"No", ElementType::No},    {"Lr", ElementType::Lr}, {"Rf", ElementType::Rf},
    {"Db", ElementType::Db},     {"Sg", ElementType::Sg},    {"Bh", ElementType::Bh}, {"Hs", ElementType::Hs},
    {"Mt", ElementType::Mt},     {"Ds", ElementType::Ds},    {"Rg", ElementType::Rg}, {"Cn", ElementType::Cn},
    {"None", ElementType::none}, {"none", ElementType::none}};

} /* namespace Utils */
} /* namespace Scine */
