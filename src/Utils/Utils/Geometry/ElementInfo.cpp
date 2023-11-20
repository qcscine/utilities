/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/ElementData.h"
#include <algorithm>

namespace Scine {
namespace Utils {
namespace detail {

// Try to separate 1H, H1, 12C, C12, C strings into alphabetical and numeric
std::pair<std::string, unsigned> isotopeInterpret(const std::string& symbol) {
  auto firstDigit = symbol.find_first_of("0123456789");
  if (firstDigit != std::string::npos) {
    // There is a digit in the string
    auto lastDigit = symbol.find_last_of("0123456789");
    unsigned atomicMass = std::stoul(symbol.substr(firstDigit, lastDigit - firstDigit + 1));

    if (lastDigit == symbol.size() - 1) {
      // Digit is postfixed
      return std::make_pair(symbol.substr(0, firstDigit), atomicMass);
    }

    // Digit is prefixed
    return std::make_pair(symbol.substr(lastDigit + 1), atomicMass);
  }

  // String contains no digits
  return std::make_pair(symbol, 0);
}

} // namespace detail

ElementType ElementInfo::elementTypeForSymbol(const std::string& symbol) {
  auto interpretPair = detail::isotopeInterpret(symbol);

  // Make all characters of the symbol lowercase if it isn't yet
  std::transform(std::begin(interpretPair.first), std::end(interpretPair.first), std::begin(interpretPair.first),
                 [](const auto c) { return std::tolower(c); });

  auto it = stringToElementType().find(interpretPair.first);
  if (it == stringToElementType().end()) {
    throw ElementSymbolNotFound(symbol);
  }

  if (interpretPair.second != 0) {
    /* There was a digit in the symbol string, so we compose an isotope using
     * the atomic mass number of the found base name and the specified atomic
     * mass number
     */
    return isotope(ElementInfo::Z(it->second), interpretPair.second);
  }

  // Return the element type as-is
  return it->second;
}

std::vector<ElementType> ElementInfo::allImplementedElements() {
  unsigned nEle = stringToElementType().size();
  std::vector<ElementType> result;
  // start at 1, because 0th is none element
  // leave out last 3, because of 2 explicit H isotope ElementTypes + electron
  for (unsigned i = 1; i < nEle - 3; ++i) {
    result.push_back(element(i));
  }
  return result;
}

std::string ElementInfo::symbol(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).symbol();
}

double ElementInfo::mass(ElementType e) {
  if (A(e) == 0) {
    return Constants::ElementDataSingleton::lookup(e).mass();
  }

  auto findIter = isotopeMap().find(e);
  if (findIter == std::end(isotopeMap())) {
    throw std::out_of_range("No data found for that isotope");
  }
  return findIter->second.mass;
}

double ElementInfo::abundance(ElementType e) {
  if (A(e) == 0) {
    throw std::logic_error("Unspecified isotope has no abundance");
  }

  auto findIter = isotopeMap().find(e);
  if (findIter == std::end(isotopeMap())) {
    throw std::out_of_range("No data found for that isotope");
  }
  return findIter->second.abundance;
}

ElementType ElementInfo::element(unsigned z) {
  return base(static_cast<ElementType>(detail::isotope(z, 0)));
}

ElementType ElementInfo::isotope(unsigned z, unsigned a) {
  auto findIter = isotopeMap().find(static_cast<ElementType>(detail::isotope(z, a)));

  if (findIter == std::end(isotopeMap())) {
    throw std::out_of_range("No such isotope!");
  }

  return findIter->first;
}

std::vector<ElementType> ElementInfo::isotopes(ElementType e) {
  std::vector<ElementType> matchingZ;
  for (const auto& iterPair : isotopeMap()) {
    if (Z(iterPair.first) == Z(e)) {
      matchingZ.push_back(iterPair.first);
    }
  }

  return matchingZ;
}

ElementType ElementInfo::base(ElementType isotope) {
  const unsigned z = Z(isotope);

  // Take care of monoisotopic elements
  switch (z) {
    case 4:
      return ElementType::Be;
    case 9:
      return ElementType::F;
    case 11:
      return ElementType::Na;
    case 13:
      return ElementType::Al;
    case 15:
      return ElementType::P;
    case 21:
      return ElementType::Sc;
    case 25:
      return ElementType::Mn;
    case 27:
      return ElementType::Co;
    case 33:
      return ElementType::As;
    case 39:
      return ElementType::Y;
    case 41:
      return ElementType::Nb;
    case 45:
      return ElementType::Rh;
    case 53:
      return ElementType::I;
    case 55:
      return ElementType::Cs;
    case 59:
      return ElementType::Pr;
    case 65:
      return ElementType::Tb;
    case 67:
      return ElementType::Ho;
    case 69:
      return ElementType::Tm;
    case 79:
      return ElementType::Au;
    case 83:
      return ElementType::Bi;
    case 87:
      return ElementType::Fr;
    case 89:
      return ElementType::Ac;
    case 91:
      return ElementType::Pa;
    case 99:
      return ElementType::Es;
    case 100:
      return ElementType::Fm;
    case 102:
      return ElementType::No;
    case 103:
      return ElementType::Lr;
    case 104:
      return ElementType::Rf;
    case 105:
      return ElementType::Db;
    case 106:
      return ElementType::Sg;
    case 107:
      return ElementType::Bh;
    case 108:
      return ElementType::Hs;
    case 109:
      return ElementType::Mt;
    case 110:
      return ElementType::Ds;
    case 111:
      return ElementType::Rg;
    case 112:
      return ElementType::Cn;
    case 113:
      return ElementType::E;
    default:
      return static_cast<ElementType>(detail::isotope(z, 0));
  }
}

double ElementInfo::paulingElectronegativity(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).paulingElectronegativity();
}

double ElementInfo::vdwRadius(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).vdWRadius();
}

double ElementInfo::covalentRadius(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).covalentRadius();
}

int ElementInfo::valElectrons(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).valElectrons();
}

int ElementInfo::sElectrons(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).sElectrons();
}

int ElementInfo::pElectrons(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).pElectrons();
}

int ElementInfo::dElectrons(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).dElectrons();
}

int ElementInfo::fElectrons(ElementType e) {
  if (A(e) != 0) {
    e = base(e);
  }

  return Constants::ElementDataSingleton::lookup(e).fElectrons();
}

const std::unordered_map<std::string, ElementType>& ElementInfo::stringToElementType() {
  static std::unordered_map<std::string, ElementType> map{
      {"none", ElementType::none}, {"h", ElementType::H},   {"d", ElementType::D},   {"t", ElementType::T},
      {"he", ElementType::He},     {"li", ElementType::Li}, {"be", ElementType::Be}, {"b", ElementType::B},
      {"c", ElementType::C},       {"n", ElementType::N},   {"o", ElementType::O},   {"f", ElementType::F},
      {"ne", ElementType::Ne},     {"na", ElementType::Na}, {"mg", ElementType::Mg}, {"al", ElementType::Al},
      {"si", ElementType::Si},     {"p", ElementType::P},   {"s", ElementType::S},   {"cl", ElementType::Cl},
      {"ar", ElementType::Ar},     {"k", ElementType::K},   {"ca", ElementType::Ca}, {"sc", ElementType::Sc},
      {"ti", ElementType::Ti},     {"v", ElementType::V},   {"cr", ElementType::Cr}, {"mn", ElementType::Mn},
      {"fe", ElementType::Fe},     {"co", ElementType::Co}, {"ni", ElementType::Ni}, {"cu", ElementType::Cu},
      {"zn", ElementType::Zn},     {"ga", ElementType::Ga}, {"ge", ElementType::Ge}, {"as", ElementType::As},
      {"se", ElementType::Se},     {"br", ElementType::Br}, {"kr", ElementType::Kr}, {"rb", ElementType::Rb},
      {"sr", ElementType::Sr},     {"y", ElementType::Y},   {"zr", ElementType::Zr}, {"nb", ElementType::Nb},
      {"mo", ElementType::Mo},     {"tc", ElementType::Tc}, {"ru", ElementType::Ru}, {"rh", ElementType::Rh},
      {"pd", ElementType::Pd},     {"ag", ElementType::Ag}, {"cd", ElementType::Cd}, {"in", ElementType::In},
      {"sn", ElementType::Sn},     {"sb", ElementType::Sb}, {"te", ElementType::Te}, {"i", ElementType::I},
      {"xe", ElementType::Xe},     {"cs", ElementType::Cs}, {"ba", ElementType::Ba}, {"la", ElementType::La},
      {"ce", ElementType::Ce},     {"pr", ElementType::Pr}, {"nd", ElementType::Nd}, {"pm", ElementType::Pm},
      {"sm", ElementType::Sm},     {"eu", ElementType::Eu}, {"gd", ElementType::Gd}, {"tb", ElementType::Tb},
      {"dy", ElementType::Dy},     {"ho", ElementType::Ho}, {"er", ElementType::Er}, {"tm", ElementType::Tm},
      {"yb", ElementType::Yb},     {"lu", ElementType::Lu}, {"hf", ElementType::Hf}, {"ta", ElementType::Ta},
      {"w", ElementType::W},       {"re", ElementType::Re}, {"os", ElementType::Os}, {"ir", ElementType::Ir},
      {"pt", ElementType::Pt},     {"au", ElementType::Au}, {"hg", ElementType::Hg}, {"tl", ElementType::Tl},
      {"pb", ElementType::Pb},     {"bi", ElementType::Bi}, {"po", ElementType::Po}, {"at", ElementType::At},
      {"rn", ElementType::Rn},     {"fr", ElementType::Fr}, {"ra", ElementType::Ra}, {"ac", ElementType::Ac},
      {"th", ElementType::Th},     {"pa", ElementType::Pa}, {"u", ElementType::U},   {"np", ElementType::Np},
      {"pu", ElementType::Pu},     {"am", ElementType::Am}, {"cm", ElementType::Cm}, {"bk", ElementType::Bk},
      {"cf", ElementType::Cf},     {"es", ElementType::Es}, {"fm", ElementType::Fm}, {"md", ElementType::Md},
      {"no", ElementType::No},     {"lr", ElementType::Lr}, {"rf", ElementType::Rf}, {"db", ElementType::Db},
      {"sg", ElementType::Sg},     {"bh", ElementType::Bh}, {"hs", ElementType::Hs}, {"mt", ElementType::Mt},
      {"ds", ElementType::Ds},     {"rg", ElementType::Rg}, {"cn", ElementType::Cn}, {"e", ElementType::E}};
  return map;
}

const std::unordered_map<ElementType, ElementInfo::IsotopeData>& ElementInfo::isotopeMap() {
  static const std::unordered_map<ElementType, ElementInfo::IsotopeData> map = {
      {ElementType::H1, {1.00782503223, 0.999885}},
      {ElementType::D, {2.01410177812, 0.000115}},
      {ElementType::T, {3.0160492779, 0.0}},
      {ElementType::He3, {3.0160293201, 0.00000134}},
      {ElementType::He4, {4.00260325413, 0.99999866}},
      {ElementType::Li6, {6.0151228874, 0.0759}},
      {ElementType::Li7, {7.0160034366, 0.9241}},
      {ElementType::Be9, {9.012183065, 1}},
      {ElementType::B10, {10.01293695, 0.199}},
      {ElementType::B11, {11.00930536, 0.801}},
      {ElementType::C12, {12.0000000, 0.9893}},
      {ElementType::C13, {13.00335483507, 0.0107}},
      {ElementType::C14, {14.0032419884, 0.0}},
      {ElementType::N14, {14.00307400443, 0.99636}},
      {ElementType::N15, {15.00010889888, 0.00364}},
      {ElementType::O16, {15.99491461957, 0.99757}},
      {ElementType::O17, {16.99913175650, 0.00038}},
      {ElementType::O18, {17.99915961286, 0.00205}},
      {ElementType::F19, {18.99840316273, 1}},
      {ElementType::Ne20, {19.9924401762, 0.9048}},
      {ElementType::Ne21, {20.993846685, 0.0027}},
      {ElementType::Ne22, {21.991385114, 0.0925}},
      {ElementType::Na23, {22.9897692820, 1}},
      {ElementType::Mg24, {23.985041697, 0.7899}},
      {ElementType::Mg25, {24.985836976, 0.1000}},
      {ElementType::Mg26, {25.982592968, 0.1101}},
      {ElementType::Al27, {26.98153853, 1}},
      {ElementType::Si28, {27.97692653465, 0.92223}},
      {ElementType::Si29, {28.97649466490, 0.04685}},
      {ElementType::Si30, {29.973770136, 0.03092}},
      {ElementType::P31, {30.97376199842, 1}},
      {ElementType::S32, {31.9720711744, 0.9499}},
      {ElementType::S33, {32.9714589098, 0.0075}},
      {ElementType::S34, {33.967867004, 0.0425}},
      {ElementType::S36, {35.96708071, 0.0001}},
      {ElementType::Cl35, {34.968852682, 0.7576}},
      {ElementType::Cl37, {36.965902602, 0.2424}},
      {ElementType::Ar36, {35.967545105, 0.003336}},
      {ElementType::Ar38, {37.96273211, 0.000629}},
      {ElementType::Ar40, {39.9623831237, 0.996035}},
      {ElementType::K39, {38.9637064864, 0.932581}},
      {ElementType::K40, {39.963998166, 0.000117}},
      {ElementType::K41, {40.9618252579, 0.067302}},
      {ElementType::Ca40, {39.962590863, 0.96941}},
      {ElementType::Ca42, {41.95861783, 0.00647}},
      {ElementType::Ca43, {42.95876644, 0.00135}},
      {ElementType::Ca44, {43.95548156, 0.02086}},
      {ElementType::Ca46, {45.9536890, 0.00004}},
      {ElementType::Ca48, {47.95252276, 0.00187}},
      {ElementType::Sc45, {44.95590828, 1}},
      {ElementType::Ti46, {45.95262772, 0.0825}},
      {ElementType::Ti47, {46.95175879, 0.0744}},
      {ElementType::Ti48, {47.94794198, 0.7372}},
      {ElementType::Ti49, {48.94786568, 0.0541}},
      {ElementType::Ti50, {49.94478689, 0.0518}},
      {ElementType::V50, {49.94715601, 0.00250}},
      {ElementType::V51, {50.94395704, 0.99750}},
      {ElementType::Cr50, {49.94604183, 0.04345}},
      {ElementType::Cr52, {51.94050623, 0.83789}},
      {ElementType::Cr53, {52.94064815, 0.09501}},
      {ElementType::Cr54, {53.93887916, 0.02365}},
      {ElementType::Mn55, {54.93804391, 1}},
      {ElementType::Fe54, {53.93960899, 0.05845}},
      {ElementType::Fe56, {55.93493633, 0.91754}},
      {ElementType::Fe57, {56.93539284, 0.02119}},
      {ElementType::Fe58, {57.93327443, 0.00282}},
      {ElementType::Co59, {58.93319429, 1}},
      {ElementType::Ni58, {57.93534241, 0.68077}},
      {ElementType::Ni60, {59.93078588, 0.26223}},
      {ElementType::Ni61, {60.93105557, 0.011399}},
      {ElementType::Ni62, {61.92834537, 0.036346}},
      {ElementType::Ni64, {63.92796682, 0.009255}},
      {ElementType::Cu63, {62.92959772, 0.6915}},
      {ElementType::Cu65, {64.92778970, 0.3085}},
      {ElementType::Zn64, {63.92914201, 0.4917}},
      {ElementType::Zn66, {65.92603381, 0.2773}},
      {ElementType::Zn67, {66.92712775, 0.0404}},
      {ElementType::Zn68, {67.92484455, 0.1845}},
      {ElementType::Zn70, {69.9253192, 0.0061}},
      {ElementType::Ga69, {68.9255735, 0.60108}},
      {ElementType::Ga71, {70.92470258, 0.39892}},
      {ElementType::Ge70, {69.92424875, 0.2057}},
      {ElementType::Ge72, {71.922075826, 0.2745}},
      {ElementType::Ge73, {72.923458956, 0.0775}},
      {ElementType::Ge74, {73.921177761, 0.3650}},
      {ElementType::Ge76, {75.921402726, 0.0773}},
      {ElementType::As75, {74.92159457, 1}},
      {ElementType::Se74, {73.922475934, 0.0089}},
      {ElementType::Se76, {75.919213704, 0.0937}},
      {ElementType::Se77, {76.919914154, 0.0763}},
      {ElementType::Se78, {77.91730928, 0.2377}},
      {ElementType::Se80, {79.9165218, 0.4961}},
      {ElementType::Se82, {81.9166995, 0.0873}},
      {ElementType::Br79, {78.9183376, 0.5069}},
      {ElementType::Br81, {80.9162897, 0.4931}},
      {ElementType::Kr78, {77.92036494, 0.00355}},
      {ElementType::Kr80, {79.91637808, 0.02286}},
      {ElementType::Kr82, {81.91348273, 0.11593}},
      {ElementType::Kr83, {82.91412716, 0.11500}},
      {ElementType::Kr84, {83.9114977282, 0.56987}},
      {ElementType::Kr86, {85.9106106269, 0.17279}},
      {ElementType::Rb85, {84.9117897379, 0.7217}},
      {ElementType::Rb87, {86.9091805310, 0.2783}},
      {ElementType::Sr84, {83.9134191, 0.0056}},
      {ElementType::Sr86, {85.9092606, 0.0986}},
      {ElementType::Sr87, {86.9088775, 0.0700}},
      {ElementType::Sr88, {87.9056125, 0.8258}},
      {ElementType::Y89, {88.9058403, 1}},
      {ElementType::Zr90, {89.9046977, 0.5145}},
      {ElementType::Zr91, {90.9056396, 0.1122}},
      {ElementType::Zr92, {91.9050347, 0.1715}},
      {ElementType::Zr94, {93.9063108, 0.1738}},
      {ElementType::Zr96, {95.9082714, 0.0280}},
      {ElementType::Nb93, {92.9063730, 1}},
      {ElementType::Mo92, {91.90680796, 0.1453}},
      {ElementType::Mo94, {93.90508490, 0.0915}},
      {ElementType::Mo95, {94.90583877, 0.1584}},
      {ElementType::Mo96, {95.90467612, 0.1667}},
      {ElementType::Mo97, {96.90601812, 0.0960}},
      {ElementType::Mo98, {97.90540482, 0.2439}},
      {ElementType::Mo100, {99.9074718, 0.0982}},
      {ElementType::Tc97, {96.9063667, 1}},
      {ElementType::Tc98, {97.9072124, 1}},
      {ElementType::Tc99, {98.9062508, 1}},
      {ElementType::Ru96, {95.90759025, 0.0554}},
      {ElementType::Ru98, {97.9052868, 0.0187}},
      {ElementType::Ru99, {98.9059341, 0.1276}},
      {ElementType::Ru100, {99.9042143, 0.1260}},
      {ElementType::Ru101, {100.9055769, 0.1706}},
      {ElementType::Ru102, {101.9043441, 0.3155}},
      {ElementType::Ru104, {103.9054275, 0.1862}},
      {ElementType::Rh103, {102.9054980, 1}},
      {ElementType::Pd102, {101.9056022, 0.0102}},
      {ElementType::Pd104, {103.9040305, 0.1114}},
      {ElementType::Pd105, {104.9050796, 0.2233}},
      {ElementType::Pd106, {105.9034804, 0.2733}},
      {ElementType::Pd108, {107.9038916, 0.2646}},
      {ElementType::Pd110, {109.90517220, 0.1172}},
      {ElementType::Ag107, {106.9050916, 0.51839}},
      {ElementType::Ag109, {108.9047553, 0.48161}},
      {ElementType::Cd106, {105.9064599, 0.0125}},
      {ElementType::Cd108, {107.9041834, 0.0089}},
      {ElementType::Cd110, {109.90300661, 0.1249}},
      {ElementType::Cd111, {110.90418287, 0.1280}},
      {ElementType::Cd112, {111.90276287, 0.2413}},
      {ElementType::Cd113, {112.90440813, 0.1222}},
      {ElementType::Cd114, {113.90336509, 0.2873}},
      {ElementType::Cd116, {115.90476315, 0.0749}},
      {ElementType::In113, {112.90406184, 0.0429}},
      {ElementType::In115, {114.903878776, 0.9571}},
      {ElementType::Sn112, {111.90482387, 0.0097}},
      {ElementType::Sn114, {113.9027827, 0.0066}},
      {ElementType::Sn115, {114.903344699, 0.0034}},
      {ElementType::Sn116, {115.90174280, 0.1454}},
      {ElementType::Sn117, {116.90295398, 0.0768}},
      {ElementType::Sn118, {117.90160657, 0.2422}},
      {ElementType::Sn119, {118.90331117, 0.0859}},
      {ElementType::Sn120, {119.90220163, 0.3258}},
      {ElementType::Sn122, {121.9034438, 0.0463}},
      {ElementType::Sn124, {123.9052766, 0.0579}},
      {ElementType::Sb121, {120.9038120, 0.5721}},
      {ElementType::Sb123, {122.9042132, 0.4279}},
      {ElementType::Te120, {119.9040593, 0.0009}},
      {ElementType::Te122, {121.9030435, 0.0255}},
      {ElementType::Te123, {122.9042698, 0.0089}},
      {ElementType::Te124, {123.9028171, 0.0474}},
      {ElementType::Te125, {124.9044299, 0.0707}},
      {ElementType::Te126, {125.9033109, 0.1884}},
      {ElementType::Te128, {127.90446128, 0.3174}},
      {ElementType::Te130, {129.906222748, 0.3408}},
      {ElementType::I127, {126.9044719, 1}},
      {ElementType::Xe124, {123.9058920, 0.000952}},
      {ElementType::Xe126, {125.9042983, 0.000890}},
      {ElementType::Xe128, {127.9035310, 0.019102}},
      {ElementType::Xe129, {128.9047808611, 0.264006}},
      {ElementType::Xe130, {129.903509349, 0.040710}},
      {ElementType::Xe131, {130.90508406, 0.212324}},
      {ElementType::Xe132, {131.9041550856, 0.269086}},
      {ElementType::Xe134, {133.90539466, 0.104357}},
      {ElementType::Xe136, {135.907214484, 0.088573}},
      {ElementType::Cs133, {132.9054519610, 1}},
      {ElementType::Ba130, {129.9063207, 0.00106}},
      {ElementType::Ba132, {131.9050611, 0.00101}},
      {ElementType::Ba134, {133.90450818, 0.02417}},
      {ElementType::Ba135, {134.90568838, 0.06592}},
      {ElementType::Ba136, {135.90457573, 0.07854}},
      {ElementType::Ba137, {136.90582714, 0.11232}},
      {ElementType::Ba138, {137.90524700, 0.71698}},
      {ElementType::La138, {137.9071149, 0.0008881}},
      {ElementType::La139, {138.9063563, 0.9991119}},
      {ElementType::Ce136, {135.90712921, 0.00185}},
      {ElementType::Ce138, {137.905991, 0.00251}},
      {ElementType::Ce140, {139.9054431, 0.88450}},
      {ElementType::Ce142, {141.9092504, 0.11114}},
      {ElementType::Pr141, {140.9076576, 1}},
      {ElementType::Nd142, {141.9077290, 0.27152}},
      {ElementType::Nd143, {142.9098200, 0.12174}},
      {ElementType::Nd144, {143.9100930, 0.23798}},
      {ElementType::Nd145, {144.9125793, 0.08293}},
      {ElementType::Nd146, {145.9131226, 0.17189}},
      {ElementType::Nd148, {147.9168993, 0.05756}},
      {ElementType::Nd150, {149.9209022, 0.05638}},
      {ElementType::Pm145, {144.9127559, 1}},
      {ElementType::Pm147, {146.9151450, 1}},
      {ElementType::Sm144, {143.9120065, 0.0307}},
      {ElementType::Sm147, {146.9149044, 0.1499}},
      {ElementType::Sm148, {147.9148292, 0.1124}},
      {ElementType::Sm149, {148.9171921, 0.1382}},
      {ElementType::Sm150, {149.9172829, 0.0738}},
      {ElementType::Sm152, {151.9197397, 0.2675}},
      {ElementType::Sm154, {153.9222169, 0.2275}},
      {ElementType::Eu151, {150.9198578, 0.4781}},
      {ElementType::Eu153, {152.9212380, 0.5219}},
      {ElementType::Gd152, {151.9197995, 0.0020}},
      {ElementType::Gd154, {153.9208741, 0.0218}},
      {ElementType::Gd155, {154.9226305, 0.1480}},
      {ElementType::Gd156, {155.9221312, 0.2047}},
      {ElementType::Gd157, {156.9239686, 0.1565}},
      {ElementType::Gd158, {157.9241123, 0.2484}},
      {ElementType::Gd160, {159.9270624, 0.2186}},
      {ElementType::Tb159, {158.9253547, 1}},
      {ElementType::Dy156, {155.9242847, 0.00056}},
      {ElementType::Dy158, {157.9244159, 0.00095}},
      {ElementType::Dy160, {159.9252046, 0.02329}},
      {ElementType::Dy161, {160.9269405, 0.18889}},
      {ElementType::Dy162, {161.9268056, 0.25475}},
      {ElementType::Dy163, {162.9287383, 0.24896}},
      {ElementType::Dy164, {163.9291819, 0.28260}},
      {ElementType::Ho165, {164.9303288, 1}},
      {ElementType::Er162, {161.9287884, 0.00139}},
      {ElementType::Er164, {163.9292088, 0.01601}},
      {ElementType::Er166, {165.9302995, 0.33503}},
      {ElementType::Er167, {166.9320546, 0.22869}},
      {ElementType::Er168, {167.9323767, 0.26978}},
      {ElementType::Er170, {169.9354702, 0.14910}},
      {ElementType::Tm169, {168.9342179, 1}},
      {ElementType::Yb168, {167.9338896, 0.00123}},
      {ElementType::Yb170, {169.9347664, 0.02982}},
      {ElementType::Yb171, {170.9363302, 0.1409}},
      {ElementType::Yb172, {171.9363859, 0.2168}},
      {ElementType::Yb173, {172.9382151, 0.16103}},
      {ElementType::Yb174, {173.9388664, 0.32026}},
      {ElementType::Yb176, {175.9425764, 0.12996}},
      {ElementType::Lu175, {174.9407752, 0.97401}},
      {ElementType::Lu176, {175.9426897, 0.02599}},
      {ElementType::Hf174, {173.9400461, 0.0016}},
      {ElementType::Hf176, {175.9414076, 0.0526}},
      {ElementType::Hf177, {176.9432277, 0.1860}},
      {ElementType::Hf178, {177.9437058, 0.2728}},
      {ElementType::Hf179, {178.9458232, 0.1362}},
      {ElementType::Hf180, {179.9465570, 0.3508}},
      {ElementType::Ta180, {179.9474648, 0.0001201}},
      {ElementType::Ta181, {180.9479958, 0.9998799}},
      {ElementType::W180, {179.9467108, 0.0012}},
      {ElementType::W182, {181.94820394, 0.2650}},
      {ElementType::W183, {182.95022275, 0.1431}},
      {ElementType::W184, {183.95093092, 0.3064}},
      {ElementType::W186, {185.9543628, 0.2843}},
      {ElementType::Re185, {184.9529545, 0.3740}},
      {ElementType::Re187, {186.9557501, 0.6260}},
      {ElementType::Os184, {183.9524885, 0.0002}},
      {ElementType::Os186, {185.9538350, 0.0159}},
      {ElementType::Os187, {186.9557474, 0.0196}},
      {ElementType::Os188, {187.9558352, 0.1324}},
      {ElementType::Os189, {188.9581442, 0.1615}},
      {ElementType::Os190, {189.9584437, 0.2626}},
      {ElementType::Os192, {191.9614770, 0.4078}},
      {ElementType::Ir191, {190.9605893, 0.373}},
      {ElementType::Ir193, {192.9629216, 0.627}},
      {ElementType::Pt190, {189.9599297, 0.00012}},
      {ElementType::Pt192, {191.9610387, 0.00782}},
      {ElementType::Pt194, {193.9626809, 0.3286}},
      {ElementType::Pt195, {194.9647917, 0.3378}},
      {ElementType::Pt196, {195.96495209, 0.2521}},
      {ElementType::Pt198, {197.9678949, 0.07356}},
      {ElementType::Au197, {196.96656879, 1}},
      {ElementType::Hg196, {195.9658326, 0.0015}},
      {ElementType::Hg198, {197.96676860, 0.0997}},
      {ElementType::Hg199, {198.96828064, 0.1687}},
      {ElementType::Hg200, {199.96832659, 0.2310}},
      {ElementType::Hg201, {200.97030284, 0.1318}},
      {ElementType::Hg202, {201.97064340, 0.2986}},
      {ElementType::Hg204, {203.97349398, 0.0687}},
      {ElementType::Tl203, {202.9723446, 0.2952}},
      {ElementType::Tl205, {204.9744278, 0.7048}},
      {ElementType::Pb204, {203.9730440, 0.014}},
      {ElementType::Pb206, {205.9744657, 0.241}},
      {ElementType::Pb207, {206.9758973, 0.221}},
      {ElementType::Pb208, {207.9766525, 0.524}},
      {ElementType::Bi209, {208.9803991, 1}},
      {ElementType::Po209, {208.9824308, 0}},
      {ElementType::Po210, {209.9828741, 0}},
      {ElementType::At210, {209.9871479, 0}},
      {ElementType::At211, {210.9874966, 0}},
      {ElementType::Rn211, {210.9906011, 0}},
      {ElementType::Rn220, {220.0113941, 0}},
      {ElementType::Rn222, {222.0175782, 0}},
      {ElementType::Fr223, {223.0197360, 1}},
      {ElementType::Ra223, {223.0185023, 0}},
      {ElementType::Ra224, {224.0202120, 0}},
      {ElementType::Ra226, {226.0254103, 0}},
      {ElementType::Ra228, {228.0310707, 0}},
      {ElementType::Ac227, {227.0277523, 1}},
      {ElementType::Th230, {230.0331341, 0}},
      {ElementType::Th232, {232.0380558, 1}},
      {ElementType::Pa231, {231.0358842, 1}},
      {ElementType::U233, {233.0396355, 0}},
      {ElementType::U234, {234.0409523, 0.000054}},
      {ElementType::U235, {235.0439301, 0.007204}},
      {ElementType::U236, {236.0455682, 0}},
      {ElementType::U238, {238.0507884, 0.992742}},
      {ElementType::Np236, {236.046570, 0}},
      {ElementType::Np237, {237.0481736, 0}},
      {ElementType::Pu238, {238.0495601, 0}},
      {ElementType::Pu239, {239.0521636, 0}},
      {ElementType::Pu240, {240.0538138, 0}},
      {ElementType::Pu241, {241.0568517, 0}},
      {ElementType::Pu242, {242.0587428, 0}},
      {ElementType::Pu244, {244.0642053, 0}},
      {ElementType::Am241, {241.0568293, 0}},
      {ElementType::Am243, {243.0613813, 0}},
      {ElementType::Cm243, {243.0613893, 0}},
      {ElementType::Cm244, {244.0627528, 0}},
      {ElementType::Cm245, {245.0654915, 0}},
      {ElementType::Cm246, {246.0672238, 0}},
      {ElementType::Cm247, {247.0703541, 0}},
      {ElementType::Cm248, {248.0723499, 0}},
      {ElementType::Bk247, {247.0703073, 0}},
      {ElementType::Bk249, {249.0749877, 0}},
      {ElementType::Cf249, {249.0748539, 0}},
      {ElementType::Cf250, {250.0764062, 0}},
      {ElementType::Cf251, {251.0795886, 0}},
      {ElementType::Cf252, {252.0816272, 0}},
      {ElementType::Es252, {252.082980, 1}},
      {ElementType::Fm257, {257.0951061, 1}},
      {ElementType::Md258, {258.0984315, 0}},
      {ElementType::Md260, {260.10365, 0}},
      {ElementType::No259, {259.10103, 1}},
      {ElementType::Lr262, {262.10961, 1}},
      {ElementType::Rf267, {267.12179, 1}},
      {ElementType::Db268, {268.12567, 1}},
      {ElementType::Sg271, {271.13393, 1}},
      {ElementType::Bh272, {272.13826, 1}},
      {ElementType::Hs270, {270.13429, 1}},
      {ElementType::Mt276, {276.15159, 1}},
      {ElementType::Ds281, {281.16451, 1}},
      {ElementType::Rg280, {280.16514, 1}},
      {ElementType::Cn285, {285.17712, 1}}};
  return map;
}

} /* namespace Utils */
} /* namespace Scine */
