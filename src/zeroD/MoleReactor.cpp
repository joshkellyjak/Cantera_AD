//! @file MoleReactor.cpp A zero-dimensional reactor with a moles as the
//! state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/utilities.h"
#include <boost/math/tools/roots.hpp>
#include <numeric>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

MoleReactor::MoleReactor(shared_ptr<Solution> sol, const string& name)
    : MoleReactor(sol, true, name)
{
}

MoleReactor::MoleReactor(shared_ptr<Solution> sol, bool clone, const string& name)
    : Reactor(sol, clone, name)
{
    m_nv = 2 + m_nsp; // internal energy, volume, and moles of each species
}

void MoleReactor::getMoles(span<CanteraDouble> y)
{
    // Use inverse molecular weights to convert to moles
    auto Y = m_thermo->massFractions();
    auto imw = m_thermo->inverseMolecularWeights();
    for (size_t i = 0; i < m_nsp; i++) {
        y[i] = m_mass * imw[i] * Y[i];
    }
}

void MoleReactor::setMassFromMoles(span<const CanteraDouble> y)
{
    auto mw = m_thermo->molecularWeights();
    // calculate mass from moles
    m_mass = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        m_mass += y[i] * mw[i];
    }
}

void MoleReactor::getState(span<CanteraDouble> y)
{
    // set the first component to the internal energy
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_thermo->intEnergy_mass() * m_mass;
    // set the second component to the total volume
    y[1] = m_vol;
    // set components y+2 ... y+K+2 to the moles of each species
    getMoles(y.subspan(m_sidx));
}

void MoleReactor::updateDefaultTolerances(span<CanteraDouble> atol, CanteraDouble baseAtol)
{
    vector<CanteraDouble> y(neq());
    getState(y);
    size_t sidx = speciesOffset();
    CanteraDouble totalMoles = std::accumulate(y.begin() + sidx,
                                        y.begin() + sidx + m_nsp, CanteraDouble(0.0));
    if (totalMoles <= 0.0) {
        return;
    }
    std::fill(atol.begin() + sidx, atol.begin() + sidx + m_nsp,
              baseAtol * totalMoles);
}

void MoleReactor::updateState(span<const CanteraDouble> y)
{
    // The components of y are [0] total internal energy, [1] the total volume, and
    // [2...K+3] are the moles of each species, and [K+3...] are the moles
    // of surface species on each wall.
    setMassFromMoles(y.subspan(m_sidx));
    m_vol = y[1];
    m_thermo->setMolesNoTruncate(y.subspan(m_sidx, m_nsp));
    if (m_energy) {
        CanteraDouble U = y[0];
        // Residual function: error in internal energy as a function of T
        auto u_err = [this, U](CanteraDouble T) {
            m_thermo->setState_TD(T, m_mass / m_vol);
            return m_thermo->intEnergy_mass() * m_mass - U;
        };
        // Boost's root finders work on plain doubles
        auto u_err_passive = [&u_err](CanteraDoublePassive T) {
            return castToPassiveDouble(u_err(T));
        };
        CanteraDouble T = m_thermo->temperature();
        boost::uintmax_t maxiter = 100;
        pair<CanteraDoublePassive, CanteraDoublePassive> TT;
        try {
            TT = bmt::bracket_and_solve_root(
                u_err_passive, castToPassiveDouble(T), 1.2, true,
                bmt::eps_tolerance<CanteraDoublePassive>(48), maxiter);
        } catch (std::exception&) {
            // Try full-range bisection if bracketing fails (for example, near
            // temperature limits for the phase's equation of state)
            try {
                TT = bmt::bisect(u_err_passive, castToPassiveDouble(m_thermo->minTemp()),
                    castToPassiveDouble(m_thermo->maxTemp()),
                    bmt::eps_tolerance<CanteraDoublePassive>(48), maxiter);
            } catch (std::exception& err2) {
                // Set m_thermo back to a reasonable state if root finding fails
                m_thermo->setState_TD(T, m_mass / m_vol);
                throw CanteraError("MoleReactor::updateState",
                    "{}\nat U = {}, rho = {}", err2.what(), U, m_mass / m_vol);
            }
        }
        if (fabs(TT.first - TT.second) > 1e-7*TT.first) {
            throw CanteraError("MoleReactor::updateState", "root finding failed");
        }
#if AD_ENABLED
        // The root was found on plain doubles. One Newton step with active values
        // gives T its derivative with respect to the state (implicit function
        // theorem: dT = -d(u_err) / (d(u_err)/dT), with u_err = 0 at the root).
        CanteraDouble resid = u_err(TT.second);
        CanteraDouble dUdT = m_thermo->cv_mass() * m_mass;
        // The corrected T has the same value as the last (passive) evaluation, so
        // drop the value-keyed property caches to recompute them with the active T
        m_thermo->invalidateCache();
        m_thermo->setState_TD(TT.second - resid / dUdT, m_mass / m_vol);
#else
        m_thermo->setState_TD(TT.second, m_mass / m_vol);
#endif
    } else {
        m_thermo->setDensity(m_mass / m_vol);
    }
    updateConnected(true);
}

void MoleReactor::eval(CanteraDouble time, span<CanteraDouble> LHS, span<CanteraDouble> RHS)
{
    auto dndt = RHS.subspan(m_sidx); // moles per time

    evalWalls(time);
    updateSurfaceProductionRates();
    // inverse molecular weights for conversion
    auto imw = m_thermo->inverseMolecularWeights();
    // volume equation
    RHS[1] = m_vdot;

    if (m_chem) {
        m_kin->getNetProductionRates(m_wdot); // "omega dot"
    }

    // Energy equation.
    // @f[
    //     \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in} - \dot m_{out} h.
    // @f]
    if (m_energy) {
        RHS[0] = - m_thermo->pressure() * m_vdot + m_Qdot;
        RHS[0] += m_thermo->intrinsicHeating() * m_vol;
    } else {
        RHS[0] = 0.0;
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dndt[k] = m_wdot[k] * m_vol + m_sdot[k];
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        // flow of species into system and dilution by other species
        for (size_t n = 0; n < m_nsp; n++) {
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
        // energy update based on mass flow
        CanteraDouble mdot = outlet->massFlowRate();
        if (m_energy) {
            RHS[0] -= mdot * m_enthalpy;
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        CanteraDouble mdot = inlet->massFlowRate();
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
        if (m_energy) {
            RHS[0] += mdot * inlet->enthalpy_mass();
        }
    }
}


size_t MoleReactor::componentIndex(const string& nm) const
{
    if (nm == "int_energy") {
        return 0;
    }
    if (nm == "volume") {
        return 1;
    }
    try {
        return m_thermo->speciesIndex(nm) + m_sidx;
    } catch (const CanteraError&) {
        throw CanteraError("MoleReactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string MoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "int_energy";
    } else if (k == 1) {
        return "volume";
    } else if (k >= m_sidx && k < neq()) {
        return m_thermo->speciesName(k - m_sidx);
    }
    throw IndexError("MoleReactor::componentName", "component", k, m_nv);
}

CanteraDouble MoleReactor::upperBound(size_t k) const {
    // Component is either int_energy, volume, or moles of a bulk or surface species
    return BigNumber;
}

CanteraDouble MoleReactor::lowerBound(size_t k) const {
    if (k == 0) {
        return -BigNumber; // int_energy
    } else if (k == 1) {
        return 0; // volume
    } else if (k >= 2 && k < m_nv) {
        return -Tiny; // moles of bulk or surface species
    } else {
        throw CanteraError("MoleReactor::lowerBound", "Index {} is out of bounds.", k);
    }
}

void MoleReactor::resetBadValues(span<CanteraDouble> y) {
    for (size_t k = m_sidx; k < m_nv; k++) {
        y[k] = std::max(y[k], 0.0);
    }
}

}
