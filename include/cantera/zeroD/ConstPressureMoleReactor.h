//! @file ConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTPRESSMOLE_REACTOR_H
#define CT_CONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/**
 * ConstPressureMoleReactor is a class for constant-pressure reactors
 * which use a state of moles.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class ConstPressureMoleReactor : public MoleReactor
{
public:
    ConstPressureMoleReactor(shared_ptr<Solution> sol, const string& name="(none)");
    ConstPressureMoleReactor(shared_ptr<Solution> sol, bool clone,
                             const string& name="(none)");

    string type() const override {
        return "ConstPressureMoleReactor";
    };

    void getState(span<CanteraDouble> y) override;
    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    void evalSteady(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    vector<size_t> initializeSteady() override;

    void updateState(span<const CanteraDouble> y) override;

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    CanteraDouble upperBound(size_t k) const override;
    CanteraDouble lowerBound(size_t k) const override;
    void resetBadValues(span<CanteraDouble> y) override;

protected:
    const size_t m_sidx = 1;
    CanteraDouble m_initialMass; //!< Initial mass [kg]; used for steady-state calculations
};

}

#endif
