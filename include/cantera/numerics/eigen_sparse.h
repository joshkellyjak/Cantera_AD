// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_SPARSE_H
#define CT_EIGEN_SPARSE_H

#include "cantera/base/config.h"
#include "cantera/base/ct_typedefs.h"
#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Sparse>
    #else
    #include <Eigen/Sparse>
    #endif
#else
#include "cantera/ext/Eigen/Sparse"
#endif

namespace Cantera
{
//! @ingroup matrices
typedef std::vector<Eigen::Triplet<CanteraDouble>> SparseTriplets;

//! @name Dense Eigen types with CanteraDouble scalars
//! Equal to Eigen's VectorXd, RowVectorXd, ArrayXd and MatrixXd in a normal build;
//! with AD enabled, their scalar type is the AD type.
//! @{
using VectorXd = Eigen::Matrix<CanteraDouble, Eigen::Dynamic, 1>;
using RowVectorXd = Eigen::Matrix<CanteraDouble, 1, Eigen::Dynamic>;
using ArrayXd = Eigen::Array<CanteraDouble, Eigen::Dynamic, 1>;
using MatrixXd = Eigen::Matrix<CanteraDouble, Eigen::Dynamic, Eigen::Dynamic>;
//! @}
}

#endif
