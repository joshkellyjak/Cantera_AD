//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * MoleReactor is meant to serve the same purpose as the reactor class but with a state
 * vector composed of moles. It also serves as the base class for other mole reactors.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class MoleReactor : public Reactor
{
public:
    MoleReactor(shared_ptr<Solution> sol, const string& name="(none)");
    MoleReactor(shared_ptr<Solution> sol, bool clone, const string& name="(none)");

    string type() const override {
        return "MoleReactor";
    }

    void getState(span<CanteraDouble> y) override;
    void updateDefaultTolerances(span<CanteraDouble> atol, CanteraDouble baseAtol) override;
    void updateState(span<const CanteraDouble> y) override;
    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override;
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    CanteraDouble upperBound(size_t k) const override;
    CanteraDouble lowerBound(size_t k) const override;
    void resetBadValues(span<CanteraDouble> y) override;

protected:
    //! Get moles of the system from mass fractions stored by thermo object
    //! @param y vector for moles to be put into
    void getMoles(span<CanteraDouble> y);

    //! Set internal mass variable based on moles given
    //! @param y vector of moles of the system
    void setMassFromMoles(span<const CanteraDouble> y);

    //! const value for the species start index
    const size_t m_sidx = 2;
};

}

#endif
