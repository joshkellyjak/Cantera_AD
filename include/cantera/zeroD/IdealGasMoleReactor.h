//! @file IdealGasMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASMOLE_REACTOR_H
#define CT_IDEALGASMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/**
 * IdealGasMoleReactor is a class for ideal gas constant-volume reactors which use a
 * state of moles.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class IdealGasMoleReactor : public MoleReactor
{
public:
    using MoleReactor::MoleReactor; // inherit constructors

    string type() const override {
        return "IdealGasMoleReactor";
    }

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    CanteraDouble upperBound(size_t k) const override;
    CanteraDouble lowerBound(size_t k) const override;
    vector<size_t> initializeSteady() override;

    void getState(span<CanteraDouble> y) override;

    void initialize(CanteraDouble t0=0.0) override;

    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    void evalSteady(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;

    void updateState(span<const CanteraDouble> y) override;

    //! Calculate an approximate Jacobian to accelerate preconditioned solvers

    //! Neglects derivatives with respect to mole fractions that would generate a
    //! fully-dense Jacobian. Connector terms are included for supported flow devices
    //! and walls, subject to derivative settings that control sparse approximations.
    void getJacobianElements(SparseTriplets& trips) override;
    void getJacobianScalingFactors(CanteraDouble& f_species, span<CanteraDouble> f_energy) override;
    void addPressureJacobian(SparseTriplets& trips, size_t row,
                             CanteraDouble coeff, bool includeSpecies=true) const override;
    void addTemperatureJacobian(SparseTriplets& trips, size_t row,
                                CanteraDouble coeff) const override;
    void addSpeciesMassFractionJacobian(SparseTriplets& trips, size_t row, size_t k,
                                        CanteraDouble coeff) const override;
    void addEnthalpyJacobian(SparseTriplets& trips, size_t row, CanteraDouble coeff,
                             bool includeComposition=true) const override;

    bool preconditionerSupported() const override {return true;};

protected:
    vector<CanteraDouble> m_uk; //!< Species molar internal energies
    CanteraDouble m_TotalCv; //!< Total heat capacity (@f$ m c_v @f$) [J/K]

     //! Initial volume [m³]; used for steady-state calculations
    CanteraDouble m_initialVolume;

     //! Initial temperature [K]; used for steady-state calculations
    CanteraDouble m_initialTemperature;
};

}

#endif
