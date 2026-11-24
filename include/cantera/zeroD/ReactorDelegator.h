//! @file ReactorDelegator.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORDELEGATOR_H
#define CT_REACTORDELEGATOR_H

#include "Reactor.h"
#include "cantera/base/Delegator.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"

namespace Cantera
{

//! An abstract base class for providing access to protected capabilities
//! Reactor objects from delegate methods, which would normally only be able to
//! access public Reactor members.
//!
//! Actual implementations of these methods are found in the templated
//! ReactorDelegator class. The purpose of this base class is so these methods
//! can be accessed by casting a Reactor* to a ReactorAccessor* without needing
//! to know the specific kind of Reactor at compilation time.
class ReactorAccessor
{
public:
    //! Set the number of equations represented by this reactor
    virtual void setNEq(size_t n) = 0;

    //! Get the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual CanteraDouble expansionRate() const = 0;

    //! Set the net rate of volume change (for example, from moving walls) [m^3/s]
    virtual void setExpansionRate(CanteraDouble v) = 0;

    //! Get the net heat transfer rate (for example, through walls) into the
    //! reactor [W]. This value is initialized and calculated as part of
    //! Reactor::evalWalls().
    virtual CanteraDouble heatRate() const = 0;

    //! Set the net heat transfer rate (for example, through walls) into the
    //! reactor [W]. For a value set using this method to affect the calculations done
    //! by Reactor::eval, this method should be called in either a "replace" or "after"
    //! delegate for Reactor::evalWalls().
    virtual void setHeatRate(CanteraDouble q) = 0;

    //! @copydoc ReactorBase::surfaceProductionRates
    virtual span<CanteraDouble> surfaceProductionRates() = 0;
};

//! Delegate methods of the Reactor class to external functions
//! @ingroup reactorGroup
template <class R>
class ReactorDelegator : public Delegator, public R, public ReactorAccessor
{
public:
    template <class... Args>
    ReactorDelegator(Args&&... args)
        : R(std::forward<Args>(args)...)
    {
        install("initialize", m_initialize, [this](CanteraDouble t0) { R::initialize(t0); });
        install("getState", m_getState,
            [this](span<CanteraDouble> y) {
                R::getState(y);
            });
        install("updateState", m_updateState,
            [this](span<CanteraDouble> y) {
                R::updateState(y);
            });
        install("updateConnected", m_updateConnected,
            [this](bool updatePressure) { R::updateConnected(updatePressure); });
        install("eval", m_eval,
            [this](CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) {
                R::eval(t, LHS, RHS);
            }
        );
        if constexpr (std::is_base_of<Reactor, R>::value) {
            install("evalWalls", m_evalWalls, [this](CanteraDouble t) { R::evalWalls(t); });
        }
        install("componentName", m_componentName,
            [this](size_t k) { return R::componentName(k); });
        install("componentIndex", m_componentIndex,
            [this](const string& nm) { return R::componentIndex(nm); });
        install("upperBound", m_upperBound,
            [this](size_t k) { return R::upperBound(k); });
        install("lowerBound", m_lowerBound,
            [this](size_t k) { return R::lowerBound(k); });
        install("getJacobianElements", m_getJacobianElements,
            [this](SparseTriplets& trips) { R::getJacobianElements(trips); });
    }

    // Overrides of Reactor methods

    string type() const override {
        return fmt::format("Extensible{}", R::type());
    }

    void initialize(CanteraDouble t0) override {
        m_initialize(t0);
    }

    void getState(span<CanteraDouble> y) override {
        checkArraySize("ReactorDelegator::getState", y.size(), R::m_nv);
        m_getState(y);
    }

    void updateState(span<const CanteraDouble> y) override {
        checkArraySize("ReactorDelegator::updateState", y.size(), R::m_nv);
        m_updateState(stripConst(y));
    }

    void updateConnected(bool updatePressure) override {
        m_updateConnected(updatePressure);
    }

    void eval(CanteraDouble t, span<CanteraDouble> LHS, span<CanteraDouble> RHS) override {
        checkArraySize("ReactorDelegator::eval[LHS]", LHS.size(), R::m_nv);
        checkArraySize("ReactorDelegator::eval[RHS]", RHS.size(), R::m_nv);
        m_eval(t, LHS, RHS);
    }

    void evalWalls(CanteraDouble t) override {
        if constexpr (std::is_base_of<Reactor, R>::value) {
            m_evalWalls(t);
        } else {
            ReactorBase::evalWalls(t);
        }
    }

    string componentName(size_t k) override {
        return m_componentName(k);
    }

    size_t componentIndex(const string& nm) const override {
        return m_componentIndex(nm);
    }

    CanteraDouble upperBound(size_t k) const override {
        return m_upperBound(k);
    }

    CanteraDouble lowerBound(size_t k) const override {
        return m_lowerBound(k);
    }

    void getJacobianElements(SparseTriplets& trips) override {
        m_getJacobianElements(trips);
    }

    // Public access to protected Reactor variables needed by derived classes

    void setNEq(size_t n) override {
        R::m_nv = n;
    }

    CanteraDouble expansionRate() const override {
        if constexpr (std::is_base_of<Reactor, R>::value) {
            return R::m_vdot;
        } else {
            throw NotImplementedError("ReactorDelegator::expansionRate",
                "Expansion rate is undefined for reactors of type '{}'.", type());
        }
    }

    void setExpansionRate(CanteraDouble v) override {
        if constexpr (std::is_base_of<Reactor, R>::value) {
            R::m_vdot = v;
        } else {
            throw NotImplementedError("ReactorDelegator::setExpansionRate",
                "Expansion rate is undefined for reactors of type '{}'.", type());
        }
    }

    CanteraDouble heatRate() const override {
        if constexpr (std::is_base_of<Reactor, R>::value) {
            return R::m_Qdot;
        } else {
            throw NotImplementedError("ReactorDelegator::heatRate",
                "Heat rate is undefined for reactors of type '{}'.", type());
        }
    }

    void setHeatRate(CanteraDouble q) override {
        if constexpr (std::is_base_of<Reactor, R>::value) {
            R::m_Qdot = q;
        } else {
            throw NotImplementedError("ReactorDelegator::setHeatRate",
                "Heat rate is undefined for reactors of type '{}'.", type());
        }
    }

    span<CanteraDouble> surfaceProductionRates() override {
        return R::m_sdot;
    }

private:
    function<void(CanteraDouble)> m_initialize;
    function<void(span<CanteraDouble>)> m_getState;
    function<void(span<CanteraDouble>)> m_updateState;
    function<void(bool)> m_updateConnected;
    function<void(CanteraDouble, span<CanteraDouble>, span<CanteraDouble>)> m_eval;
    function<void(CanteraDouble)> m_evalWalls;
    function<string(size_t)> m_componentName;
    function<size_t(const string&)> m_componentIndex;
    function<CanteraDouble(size_t)> m_upperBound;
    function<CanteraDouble(size_t)> m_lowerBound;
    function<void(SparseTriplets&)> m_getJacobianElements;
};

}
#endif
