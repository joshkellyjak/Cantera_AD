//! @file IdealGasConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/ConstPressureMoleReactor.h"

namespace Cantera
{

/**
 * IdealGasConstPressureMoleReactor is a class for constant-pressure reactors
 * which use temperature and species moles as state variables. The class name is
 * historical; this formulation is applicable to non-ideal equations of state
 * where the ThermoPhase implements the required enthalpy and heat capacity
 * properties.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class IdealGasConstPressureMoleReactor : public ConstPressureMoleReactor
{
public:
    using ConstPressureMoleReactor::ConstPressureMoleReactor; // inherit constructors

    string type() const override {
        return "IdealGasConstPressureMoleReactor";
    };

    void getState(span<CanteraDouble> y) override;

    void initialize(CanteraDouble t0=0.0) override;

    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;

    void updateState(span<const CanteraDouble> y) override;
    void getJacobianScalingFactors(CanteraDouble& f_species, span<CanteraDouble> f_energy) override;

    //! Calculate an approximate Jacobian to accelerate preconditioned solvers
    //!
    //! Neglects derivatives with respect to mole fractions that would generate a
    //! fully-dense Jacobian. Connector terms are included for supported flow devices
    //! and walls, subject to derivative settings that control sparse approximations.
    void getJacobianElements(SparseTriplets& trips) override;
    void addTemperatureJacobian(SparseTriplets& trips, size_t row,
                                CanteraDouble coeff) const override;
    void addSpeciesMassFractionJacobian(SparseTriplets& trips, size_t row, size_t k,
                                        CanteraDouble coeff) const override;
    void addEnthalpyJacobian(SparseTriplets& trips, size_t row, CanteraDouble coeff,
                             bool includeComposition=true) const override;

    bool preconditionerSupported() const override { return true; };

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    CanteraDouble upperBound(size_t k) const override;
    CanteraDouble lowerBound(size_t k) const override;

protected:
    vector<CanteraDouble> m_hk; //!< Species molar enthalpies
    CanteraDouble m_TotalCp; //!< Total heat capacity (@f$ m c_p @f$) [J/K]
};

}

#endif
