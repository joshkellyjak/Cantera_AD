//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTP_REACTOR_H
#define CT_IDEALGASCONSTP_REACTOR_H

#include "ConstPressureReactor.h"

namespace Cantera
{

/**
 * Class IdealGasConstPressureReactor is a class for constant-pressure reactors
 * that use temperature instead of total enthalpy as a state variable. The class
 * name is historical; this formulation is applicable to non-ideal equations of
 * state where the ThermoPhase implements the required enthalpy and heat
 * capacity properties.
 * @ingroup reactorGroup
 */
class IdealGasConstPressureReactor : public ConstPressureReactor
{
public:
    using ConstPressureReactor::ConstPressureReactor; // inherit constructors

    string type() const override {
        return "IdealGasConstPressureReactor";
    }

    void getState(span<CanteraDouble> y) override;

    void initialize(CanteraDouble t0=0.0) override;
    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    void evalSteady(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    void updateState(span<const CanteraDouble> y) override;
    vector<size_t> initializeSteady() override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "temperature", the name of a homogeneous phase species, or the name of a
    //! surface species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    CanteraDouble upperBound(size_t k) const override;
    CanteraDouble lowerBound(size_t k) const override;

protected:
    vector<CanteraDouble> m_hk; //!< Species molar enthalpies

     //! Initial mass [kg]; used for steady-state calculations
    CanteraDouble m_initialMass;

    //! Initial temperature [K]; used for steady-state calculations
    CanteraDouble m_initialTemperature;
};
}

#endif
