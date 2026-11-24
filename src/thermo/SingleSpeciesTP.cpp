/**
 *  @file SingleSpeciesTP.cpp
 *  Definitions for the SingleSpeciesTP class, which is a filter class for ThermoPhase,
 *  that eases the construction of single species phases
 *  ( see @ref thermoprops and class @link Cantera::SingleSpeciesTP SingleSpeciesTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SingleSpeciesTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{

// ------------ Molar Thermodynamic Properties --------------------

CanteraDouble SingleSpeciesTP::enthalpy_mole() const
{
    CanteraDouble hbar;
    getPartialMolarEnthalpies(span<CanteraDouble>(&hbar, 1));
    return hbar;
}

CanteraDouble SingleSpeciesTP::intEnergy_mole() const
{
    CanteraDouble ubar;
    getPartialMolarIntEnergies(span<CanteraDouble>(&ubar, 1));
    return ubar;
}

CanteraDouble SingleSpeciesTP::entropy_mole() const
{
    CanteraDouble sbar;
    getPartialMolarEntropies(span<CanteraDouble>(&sbar, 1));
    return sbar;
}

CanteraDouble SingleSpeciesTP::gibbs_mole() const
{
    CanteraDouble gbar;

    // Get the chemical potential of the first species. This is the same as the
    // partial molar Gibbs free energy.
    getChemPotentials(span<CanteraDouble>(&gbar, 1));
    return gbar;
}

CanteraDouble SingleSpeciesTP::cp_mole() const
{
    CanteraDouble cpbar;

    // Really should have a partial molar heat capacity function in ThermoPhase.
    // However, the standard state heat capacity will do fine here for now.
    getCp_R(span<CanteraDouble>(&cpbar, 1));
    cpbar *= GasConstant;
    return cpbar;
}

CanteraDouble SingleSpeciesTP::cv_mole() const
{
    // For single species, we go directory to the general Cp - Cv relation
    //
    //     Cp = Cv + alpha**2 * V * T / beta
    //
    // where
    //     alpha = volume thermal expansion coefficient
    //     beta  = isothermal compressibility
    CanteraDouble cvbar = cp_mole();
    CanteraDouble alpha = thermalExpansionCoeff();
    CanteraDouble beta = isothermalCompressibility();
    CanteraDouble V = molecularWeight(0)/density();
    CanteraDouble T = temperature();
    if (beta != 0.0) {
        cvbar -= alpha * alpha * V * T / beta;
    }
    return cvbar;
}

// ----------- Partial Molar Properties of the Solution -----------------

void SingleSpeciesTP::getChemPotentials(span<CanteraDouble> mu) const
{
    getStandardChemPotentials(mu);
}

void SingleSpeciesTP::getPartialMolarEnthalpies(span<CanteraDouble> hbar) const
{
    getEnthalpy_RT(hbar);
    hbar[0] *= RT();
}

void SingleSpeciesTP::getPartialMolarIntEnergies(span<CanteraDouble> ubar) const
{
    getIntEnergy_RT(ubar);
    ubar[0] *= RT();
}

void SingleSpeciesTP::getPartialMolarEntropies(span<CanteraDouble> sbar) const
{
    getEntropy_R(sbar);
    sbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarCp(span<CanteraDouble> cpbar) const
{
    getCp_R(cpbar);
    cpbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarVolumes(span<CanteraDouble> vbar) const
{
    checkArraySize("SingleSpeciesTP::getPartialMolarVolumes", vbar.size(), 1);
    vbar[0] = molecularWeight(0) / density();
}

// Properties of the Standard State of the Species in the Solution

void SingleSpeciesTP::getStandardVolumes(span<CanteraDouble> vbar) const
{
    checkArraySize("SingleSpeciesTP::getStandardVolumes", vbar.size(), 1);
    vbar[0] = molecularWeight(0) / density();
}

// ---- Thermodynamic Values for the Species Reference States -------

void SingleSpeciesTP::getEnthalpy_RT_ref(span<CanteraDouble> hrt) const
{
    checkArraySize("SingleSpeciesTP::getEnthalpy_RT_ref", hrt.size(), 1);
    _updateThermo();
    hrt[0] = m_h0_RT;
}

void SingleSpeciesTP::getGibbs_RT_ref(span<CanteraDouble> grt) const
{
    checkArraySize("SingleSpeciesTP::getGibbs_RT_ref", grt.size(), 1);
    _updateThermo();
    grt[0] = m_h0_RT - m_s0_R;
}

void SingleSpeciesTP::getGibbs_ref(span<CanteraDouble> g) const
{
    getGibbs_RT_ref(g);
    g[0] *= RT();
}

void SingleSpeciesTP::getEntropy_R_ref(span<CanteraDouble> er) const
{
    checkArraySize("SingleSpeciesTP::getEntropy_R_ref", er.size(), 1);
    _updateThermo();
    er[0] = m_s0_R;
}

void SingleSpeciesTP::getCp_R_ref(span<CanteraDouble> cpr) const
{
    checkArraySize("SingleSpeciesTP::getCp_R_ref", cpr.size(), 1);
    _updateThermo();
    cpr[0] = m_cp0_R;
}

bool SingleSpeciesTP::addSpecies(shared_ptr<Species> spec)
{
    if (m_kk != 0) {
        throw CanteraError("SingleSpeciesTP::addSpecies",
            "Stoichiometric substances may only contain one species.");
    }
    return ThermoPhase::addSpecies(spec);
}

void SingleSpeciesTP::_updateThermo() const
{
    CanteraDouble tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo.update_single(0, tnow, m_cp0_R, m_h0_RT, m_s0_R);
        m_tlast = tnow;
    }
}

}
