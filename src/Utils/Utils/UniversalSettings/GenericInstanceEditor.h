/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef UNIVERSALSETTINGS_GENERICINSTANCEEDITOR_H
#define UNIVERSALSETTINGS_GENERICINSTANCEEDITOR_H
/* Internal Headers */
#include "ValueCollection.h"
/* External Headers */
#include <cassert>
#include <memory>

namespace Scine {
namespace Utils {
namespace UniversalSettings {
/* Forward Declarations */
class DescriptorCollection;

/*!
 * Same as GenericInstanceEditor, defined below, without create function for the case that the constructor
 * needs more parameters, which can then be specified in the derived class.
 */
template<typename BaseT>
class GenericInstanceEditor {
 public:
  using BaseClass = BaseT;

  virtual ~GenericInstanceEditor() = default;

  /*! Get the descriptors for the settings for BaseClass. */
  virtual DescriptorCollection getSettingDescriptors() const = 0;

  /*! Get the current parameters for the given instance. */
  virtual ValueCollection getAppliedSettings(const BaseClass& instance) const = 0;

  /*! Enquire whether the Editor class can handle the given instance. */
  virtual bool relatesToInstance(const BaseClass& instance) const = 0;

  /*! Apply settings to some BaseClass instance. */
  virtual void apply(BaseClass& instance, const ValueCollection& values) const = 0;
};

/*!
 * Template for a class being able to create and modify instances of some polymorphic type from settings
 * specified in the UniversalSettings syntax.
 */
template<typename Base>
class GenericInstanceEditorWithDefaultConstructor : public GenericInstanceEditor<Base> {
 public:
  using BaseClass = Base;

  /*! Create an instance with given setting values. */
  virtual std::unique_ptr<Base> create(const ValueCollection& values) const = 0;

  /*! Create an instance with default values. */
  std::unique_ptr<Base> createDefault() const {
    auto defaultValues = createDefaultValueCollection(this->getSettingDescriptors());
    return create(defaultValues);
  }
};

/*!
 * Specification of GenericInstanceEditorWithoutCreateFunction with the type T, to generate automatically some of the
 * virtual functions.
 * \tparam BaseEditor must be some GenericInstanceEditor
 * \tparam T class to instantiate, must be a * derived class of the template parameter of BaseEditor.
 */
template<typename BaseEditor, typename T>
class GenericInstanceEditorImpl : public BaseEditor {
 public:
  using InstanceClass = T;
  using BaseClass = typename BaseEditor::BaseClass;

  bool relatesToInstance(const BaseClass& instance) const final {
    const auto* derivedPtr = dynamic_cast<T const*>(&instance);
    return derivedPtr != nullptr;
  }

  void apply(BaseClass& instance, const ValueCollection& values) const final {
    auto& derived = castDown(instance);
    applyImpl(derived, values);
  }

  ValueCollection getAppliedSettings(const BaseClass& instance) const final {
    const auto& derived = castDown(instance);
    return getAppliedSettingsImpl(derived);
  }

 protected:
  T& castDown(BaseClass& instance) const {
    assert(relatesToInstance(instance));
    auto& derivedInstance = dynamic_cast<T&>(instance); // throws std::bad_cast if not possible
    return derivedInstance;
  }
  const T& castDown(const BaseClass& instance) const {
    assert(relatesToInstance(instance));
    const auto& derivedInstance = dynamic_cast<const T&>(instance); // throws std::bad_cast if not possible
    return derivedInstance;
  }

 private:
  virtual void applyImpl(T& instance, const ValueCollection& values) const = 0;
  virtual ValueCollection getAppliedSettingsImpl(const T& instance) const = 0;
};

/*!
 * Specification of GenericInstanceEditor with the type T, to generate automatically some of the virtual functions.
 * \tparam BaseEditor must be some GenericInstanceEditorWithDefaultConstructor
 * \tparam T class to instantiate, must be a derived class of the template parameter of BaseEditor.
 */
template<typename BaseEditor, typename T>
class GenericInstanceEditorWithDefaultConstructorImpl : public GenericInstanceEditorImpl<BaseEditor, T> {
 public:
  using InstanceClass = T;
  using BaseClass = typename BaseEditor::BaseClass;

  std::unique_ptr<BaseClass> create(const ValueCollection& values) const override {
    auto instance = std::make_unique<T>();
    this->apply(*instance, values);
    return instance;
  }
};

} /* namespace UniversalSettings */
} /* namespace Utils */
} /* namespace Scine */

#endif // UNIVERSALSETTINGS_GENERICINSTANCEEDITOR_H
