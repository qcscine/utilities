/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef SCINE_DFTD3PARAMETERS_H
#define SCINE_DFTD3PARAMETERS_H

#include <array>
#include <vector>

namespace Scine {
namespace Utils {

enum class ElementType : unsigned;

namespace Dftd3 {

/**
 * @class Dftd3Parameters Dftd3Parameters.h
 * @brief A class handling the parameters involved in the calculation of D3 energies and gradients.
 *
 *      Most parameters in this method are static constexpr because they stay the same for all D3 methods.
 *      The parameters a1, s8 and a2 are specified by the user.
 *      One set of parameters, the reference pairs used to calculate the C6 coefficients, are outsourced to the
 *      Dftd3ReferencePairs class, because it involves a very large number of values.
 */
class Dftd3Parameters {
 public:
  /**
   * @brief Default constructor.
   */
  Dftd3Parameters() = default;
  /**
   * @brief Constructor that allows for passing a1, s8 and a2 parameters.
   */
  Dftd3Parameters(double a1, double s8, double a2);
  /**
   * @brief Default destructor.
   */
  ~Dftd3Parameters() = default;

  /**
   * @brief Getter for a covalent radius for a given element type.
   * @return double
   */
  double getCovalentRadius(ElementType elementType);
  /**
   * @brief Getter for the k1 parameter.
   * @return double
   */
  double getK1();
  /**
   * @brief Getter for the k3 parameter.
   * @return double
   */
  double getK3();
  /**
   * @brief Getter for reference pairs for one atom pair.
   * @return const std::array<std::array<double, 3>, 75>& An 75 x 3 array which contains all the reference pairs for the
   * atom pair consisting of element type 1 and element type 2. Some of the reference pairs are just placeholders.
   */
  const std::array<std::array<double, 3>, 25>& getReferencePairs(ElementType elementType1, ElementType elementType2);

  /**
   * @brief Getter for the r2r4 parameter for a given element type.
   *
   *     The name r2r4 is copied from Grimme's implementation of D3.
   *     It refers to the square root of the value Q, with Q being defined in
   *     Eq.(9) of J. Chem. Phys. 132, 154104, (2010).
   *
   * @return double
   */
  double getR2r4(ElementType elementType);
  /**
   * @brief Getter for the a1 parameter.
   * @return double
   */
  double getA1();
  /**
   * @brief Getter for the s8 parameter.
   * @return double
   */
  double getS8();
  /**
   * @brief Getter for the a1 parameter.
   * @return double
   */
  double getA2();

 private:
  // These three parameters have to be specified by the user as they depend on the method chosen
  double a1_;
  double s8_;
  double a2_;

  // All the following parameters for have been taken from Grimme's D3 implementation
  static constexpr double covalentRadii_[94] = {
      0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, 1.88972601, 1.78894056, 1.58736983, 1.61256616,
      1.68815527, 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, 2.57002732, 2.49443835, 2.41884923,
      4.43455700, 3.88023730, 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, 2.62041997, 2.51963467,
      2.49443835, 2.54483100, 2.74640188, 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, 2.94797246,
      4.76210950, 4.20778980, 3.70386304, 3.50229216, 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717,
      2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, 3.17473967, 3.09915070, 3.32591790, 3.30072128,
      5.26603625, 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, 3.93062995, 3.90543362, 3.80464833,
      3.82984466, 3.80464833, 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, 3.67866672, 3.45189952,
      3.30072128, 3.09915070, 2.97316878, 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, 3.27552496,
      3.27552496, 3.42670319, 3.30072128, 3.47709584, 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289,
      3.82984466, 3.85504098, 3.88023730, 3.90543362};

  static constexpr double r2r4_[94] = {
      2.00734898,  1.56637132,  5.01986934, 3.85379032, 3.64446594, 3.10492822,  2.71175247,  2.59361680, 2.38825250,
      2.21522516,  6.58585536,  5.46295967, 5.65216669, 4.88284902, 4.29727576,  4.04108902,  3.72932356, 3.44677275,
      7.97762753,  7.07623947,  6.60844053, 6.28791364, 6.07728703, 5.54643096,  5.80491167,  5.58415602, 5.41374528,
      5.28497229,  5.22592821,  5.09817141, 6.12149689, 5.54083734, 5.06696878,  4.87005108,  4.59089647, 4.31176304,
      9.55461698,  8.67396077,  7.97210197, 7.43439917, 6.58711862, 6.19536215,  6.01517290,  5.81623410, 5.65710424,
      5.52640661,  5.44263305,  5.58285373, 7.02081898, 6.46815523, 5.98089120,  5.81686657,  5.53321815, 5.25477007,
      11.02204549, 10.15679528, 9.35167836, 9.06926079, 8.97241155, 8.90092807,  8.85984840,  8.81736827, 8.79317710,
      7.89969626,  8.80588454,  8.42439218, 8.54289262, 8.47583370, 8.45090888,  8.47339339,  7.83525634, 8.20702843,
      7.70559063,  7.32755997,  7.03887381, 6.68978720, 6.05450052, 5.88752022,  5.70661499,  5.78450695, 7.79780729,
      7.26443867,  6.78151984,  6.67883169, 6.39024318, 6.09527958, 11.79156076, 11.10997644, 9.51377795, 8.67197068,
      8.77140725,  8.65402716,  8.53923501, 8.85024712};

  static constexpr double k1_ = 16.0;
  static constexpr double k3_ = 4.0;
};

} // namespace Dftd3
} // namespace Utils
} // namespace Scine

#endif // SCINE_DFTD3PARAMETERS_H
