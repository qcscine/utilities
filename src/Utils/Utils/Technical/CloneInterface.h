/**
 * @file CloneInterface.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UTILS_CLONEINTERFACE_H
#define UTILS_CLONEINTERFACE_H

#include <memory>

namespace Scine {
namespace Utils {

/**
 * @brief Strong typing for the abstract level class.
 * To provide template parameter resolution at compile time, a strong type
 * is provided to differentiate the leaf class from the abstract, middle
 * class.
 */
template<class T>
class Abstract {};

/**
 * @brief Class serving the purpose of a generic clone interface.
 * This class provides the clone method for any interface in core needing it.
 * It uses the curiously recurrent template pattern (CRTP) to acquire infos
 * at compile time about both the interface and the derived class.
 * This allows for the one-time-only generation of code for the clone method,
 * avoiding unnecessary code duplication.
 * Basically, a derived class must just inherit from the CloneInterface with
 * template parameters Derived = derived class and Base = interface class to
 * have a functioning clone method.
 * See
 * https://www.fluentcpp.com/2017/09/12/how-to-return-a-smart-pointer-and-use-covariance/
 * for further details.
 *
 * Note that classes deriving from CloneInterface<Derived, Base> do implicitly
 * derive from Base, hence they do not need an additional 'public Base'
 * statement.
 */
template<class Derived, class Base, class Interface = Base>
class CloneInterface : public Base {
 public:
  ~CloneInterface() override = default;
  /**
   * @brief Templetized clone method, it hides the interface clone method.
   * @tparam Derived The type of the derived class.
   * @return A unique_ptr to the derived class.
   * This function allows for working directly with the derived class if the
   * type is known. It hides the Base class' clone() method (does not override
   * it!).
   * This is achieved with templetized Derived and Base classes.
   */
  std::shared_ptr<Derived> clone() const {
    return std::dynamic_pointer_cast<Derived>(this->cloneImpl());
  }

 private:
  /*
   * This cloneImpl() overrides the Base class' method allowing for
   * covariant return type in the base class' clone().
   */
  std::shared_ptr<Interface> cloneImpl() const override {
    return std::make_shared<Derived>((static_cast<const Derived&>(*this)));
  }
};

/**
 * @brief Specialization for abstract classes.
 * This is necessary, because otherwise the inheritance pattern would not work.
 * It would not be necessary if the inheritance would just have 2 level:
 * interface class and derived class. From the moment that there are also
 * abstract classes between interface and derived class, as the Core::Calculator
 * in the Sparrow module is the most notorious example, the cloneImpl() method
 * must be virtualized in the abstract class.
 */
template<class Derived, class Base, class Interface>
class CloneInterface<Abstract<Derived>, Base, Interface> : public Base {
 public:
  ~CloneInterface() override = default;

  std::shared_ptr<Derived> clone() const {
    return std::dynamic_pointer_cast<Derived>(this->cloneImpl());
  }

 private:
  std::shared_ptr<Interface> cloneImpl() const override = 0;
};

} // namespace Utils
} // namespace Scine

#endif // UTILS_CLONEINTERFACE_H
