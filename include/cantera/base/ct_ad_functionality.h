#pragma once
#include "ct_typedefs.h"
#include "fmt.h"

#if AD_ENABLED

#include <cmath>

//! Strip the derivative information from AD values; pass all other arguments
//! (strings, integers, ...) through unchanged. Used for printf-style formatting,
//! which cannot handle CoDiPack types.
template <typename T>
inline decltype(auto) adPassiveArg(T const& v) {
    if constexpr (codi::ExpressionTraits::IsExpression<T>::value) {
        return codi::RealTraits::getPassiveValue(v);
    } else {
        return (v);
    }
}

template <typename ... Args>
inline std::string sprintfOverload(char const* message, Args&& ... args) {

    return fmt::sprintf(message, adPassiveArg(args)...);
}

//! expm1 for AD types (CoDiPack does not provide it). The primal value is computed
//! with std::expm1 to keep its accuracy for small arguments; the derivative is
//! d/dx expm1(x) = exp(x).
template <typename Expr, typename = codi::ExpressionTraits::EnableIfExpression<Expr>>
inline CanteraDouble expm1(Expr const& arg) {
    const CanteraDouble x = arg;
    const double v = x.getValue();
    if (v > 700.0) {
        return exp(x) - 1.0;
    }
    return std::expm1(v) + std::exp(v) * (x - v);
}

template<typename To, typename From, typename = void>
struct ADStaticCastImpl {
    inline static To cast(From const& from) {
        return static_cast<To>(from);
    }
};

template<typename To, typename From>
struct ADStaticCastImpl<To, From, codi::ExpressionTraits::EnableIfExpression<From>> {
    inline static To cast(From const& from) {
        return static_cast<To>(from.getValue());
    }
};

template<typename To, typename From>
inline To ad_static_cast(From const& from) {
    return ADStaticCastImpl<To, From>::cast(from);
}
#define cantera_cast ad_static_cast

template <> struct fmt::formatter<CanteraDouble> : fmt::formatter<double> {

  auto format(const CanteraDouble& p, format_context& ctx) const {
    return fmt::formatter<double>::format(p.getValue(), ctx);
  }
};

template <typename R, template<typename> class Op, typename ... Exps> struct fmt::formatter<codi::ComputeExpression<R, Op, Exps...>> : fmt::formatter<double> {

  auto format(const CanteraDouble& p, format_context& ctx) const {
    return fmt::formatter<double>::format(p.getValue(), ctx);
  }
};

inline CanteraDoublePassive castToPassiveDouble(CanteraDouble const& value) {
  return value.getValue();
}

inline std::vector<CanteraDouble> castToDoubleVector(std::vector<CanteraDoublePassive> const& vec) {
  std::vector<CanteraDouble> res(vec.size());
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = vec[i];
  }
  return res;
}

inline std::vector<std::vector<CanteraDouble>> castToDoubleVectorVector(std::vector<std::vector<CanteraDoublePassive>> const& vec) {
  std::vector<std::vector<CanteraDouble>> res(vec.size(), std::vector<CanteraDouble>(0));
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = castToDoubleVector(vec[i]);
  }
  return res;
}

inline std::vector<CanteraDoublePassive> castToPassiveDoubleVector(std::vector<CanteraDouble> const& vec) {
  std::vector<CanteraDoublePassive> res(vec.size());
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = vec[i].getValue();
  }
  return res;
}

inline std::vector<std::vector<CanteraDoublePassive>> castToPassiveDoubleVectorVector(std::vector<std::vector<CanteraDouble>> const& vec) {
  std::vector<std::vector<CanteraDoublePassive>> res(vec.size(), std::vector<CanteraDoublePassive>(0));
  
  for(size_t i = 0; i < vec.size(); i += 1) {
    res[i] = castToPassiveDoubleVector(vec[i]);
  }
  return res;
}
#else
#define cantera_cast static_cast
#define sprintfOverload fmt::sprintf

inline CanteraDoublePassive castToPassiveDouble(CanteraDouble const& value) { return value; }

inline std::vector<CanteraDouble> const& castToDoubleVector(std::vector<CanteraDoublePassive> const& vec) {
  return vec;
}

inline std::vector<std::vector<CanteraDouble>> const& castToDoubleVectorVector(std::vector<std::vector<CanteraDoublePassive>> const& vec) {
  return vec;
}
inline std::vector<CanteraDoublePassive> const& castToPassiveDoubleVector(std::vector<CanteraDouble> const& vec) {
  return vec;
}
inline std::vector<std::vector<CanteraDoublePassive>> const& castToPassiveDoubleVectorVector(std::vector<std::vector<CanteraDouble>> const& vec) {
  return vec;
}
#endif
