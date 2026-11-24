// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_DENSE_H
#define CT_EIGEN_DENSE_H

#include "cantera/base/ct_defs.h"

#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Dense>
    #else
    #include <Eigen/Dense>
    #endif
#else
#include "cantera/ext/Eigen/Dense"
#endif
#include <concepts>

namespace Cantera
{

//! @addtogroup matrices
//! @{

//! @name Dense Eigen types with CanteraDouble scalars
//! Equal to Eigen's VectorXd, RowVectorXd, ArrayXd and MatrixXd in a normal build;
//! with AD enabled, their scalar type is the AD type.
//! @{
using VectorXd = Eigen::Matrix<CanteraDouble, Eigen::Dynamic, 1>;
using RowVectorXd = Eigen::Matrix<CanteraDouble, 1, Eigen::Dynamic>;
using ArrayXd = Eigen::Array<CanteraDouble, Eigen::Dynamic, 1>;
using MatrixXd = Eigen::Matrix<CanteraDouble, Eigen::Dynamic, Eigen::Dynamic>;
//! @}

typedef Eigen::Map<MatrixXd> MappedMatrix;
typedef Eigen::Map<const MatrixXd> ConstMappedMatrix;
typedef Eigen::Map<VectorXd> MappedVector;
typedef Eigen::Map<const VectorXd> ConstMappedVector;
typedef Eigen::Map<RowVectorXd> MappedRowVector;
typedef Eigen::Map<const RowVectorXd> ConstMappedRowVector;

//! @}

template<class Derived>
concept EigenDenseDouble = std::same_as<typename Derived::Scalar, CanteraDouble>;

// Require a 1D type with contiguous storage
template<class Derived>
concept EigenVectorLike = EigenDenseDouble<Derived>
    && (Derived::IsVectorAtCompileTime == 1);

//! Convenience wrapper for accessing Eigen vector/array/map data as a span
//! @todo Remove once %Cantera requires Eigen 5.0 or newer.
template<EigenVectorLike Derived>
inline span<CanteraDouble> asSpan(Eigen::DenseBase<Derived>& v)
{
    return span<CanteraDouble>(v.derived().data(), v.size());
}

//! Convenience wrapper for accessing Eigen vector/array/map data as a span
//! @todo Remove once %Cantera requires Eigen 5.0 or newer.
template<EigenVectorLike Derived>
inline span<const CanteraDouble> asSpan(const Eigen::DenseBase<Derived>& v)
{
    return span<const CanteraDouble>(v.derived().data(), v.size());
}

//! Convenience wrapper for accessing std::vector as an Eigen VectorXd
inline MappedVector asVectorXd(vector<CanteraDouble>& v)
{
    return MappedVector(v.data(), static_cast<Eigen::Index>(v.size()));
}

//! Convenience wrapper for accessing std::vector as an Eigen VectorXd
inline ConstMappedVector asVectorXd(const vector<CanteraDouble>& v)
{
    return ConstMappedVector(v.data(), static_cast<Eigen::Index>(v.size()));
}

//! Convenience wrapper for accessing std::span as an Eigen VectorXd
inline MappedVector asVectorXd(span<CanteraDouble> v)
{
    return MappedVector(v.data(), static_cast<Eigen::Index>(v.size()));
}

//! Convenience wrapper for accessing std::span as an Eigen VectorXd
inline ConstMappedVector asVectorXd(span<const CanteraDouble> v)
{
    return ConstMappedVector(v.data(), static_cast<Eigen::Index>(v.size()));
}

}
#endif
