#ifndef INITIALIZER_INPUTAUX_H_
#define INITIALIZER_INPUTAUX_H_

#include <type_traits>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace initializers {

    /**
     * Returns a value from Yaml::Node converting it to from a string to the type derived from DefaultValue
     * */
    template <typename T>
    T getParamIfExists(const YAML::Node& Param, std::string&& Field, T DefaultValue) {
      if (std::is_same<T, bool>::value) {
        T Value{DefaultValue};
        if (Param[Field]) {
          Value = Param[Field].as<int>() > 0;
        }
        return Value;
      }
      else {
        return Param[Field] ? Param[Field].as<T>() : DefaultValue;
      }
    }

    /**
     * \brief Returns true if number elements in the input string (separated by the white space)
     *  is less or equal to the size of a container
     * */
    template<typename OutputType, typename ContainerT>
    bool isCapacityEnough(const std::string& InputString, ContainerT& OutputMask) {
      std::istringstream InputStream(InputString);
      auto Begin = std::istream_iterator<OutputType>(InputStream);
      auto End = std::istream_iterator<OutputType>();

      const int NumInputElements = std::distance(Begin, End);
      return NumInputElements <= OutputMask.size();
    }

    /**
     * \brief Initializes the given input mask with values from the input string
     *
     * \throws runtime_error if an input string contains more parameters than the capacity of a provided container
     * */
    template<typename ContainerT>
    void convertStringToMask(const std::string& StringMask, ContainerT& Mask) {
      using T = typename std::iterator_traits<typename ContainerT::iterator>::value_type;

      if (!isCapacityEnough<T>(StringMask, Mask))
        throw std::runtime_error("Num. input elements is more than the Mask capacity");

      std::istringstream InputStream(StringMask);
      auto Begin = std::istream_iterator<T>(InputStream);
      auto End = std::istream_iterator<T>();

      for (int Index = 0; Begin != End; ++Index, ++Begin) {
        if (std::is_same<T, bool>::value) {
          Mask[Index] = (*Begin) > 0;
        }
        else {
          Mask[Index] = (*Begin);
        }
      }
    }

  }
}

#endif //INITIALIZER_INPUTAUX_H_
