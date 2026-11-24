/**
 *  @file IdealSolnGasVPSS.cpp
 * Definition file for a derived class of ThermoPhase that assumes
 * an ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see @ref thermoprops and
 * class @link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

IdealSolnGasVPSS::IdealSolnGasVPSS(const string& infile, string id_)
{
    initThermoFile(infile, id_);
}

void IdealSolnGasVPSS::setStandardConcentrationModel(const string& model)
{
    if (caseInsensitiveEquals(model, "unity")) {
        m_formGC = 0;
    } else if (caseInsensitiveEquals(model, "species-molar-volume")
               || caseInsensitiveEquals(model, "molar_volume")) {
        m_formGC = 1;
    } else if (caseInsensitiveEquals(model, "solvent-molar-volume")
               || caseInsensitiveEquals(model, "solvent_volume")) {
        m_formGC = 2;
    } else {
        throw CanteraError("IdealSolnGasVPSS::setStandardConcentrationModel",
                           "Unknown standard concentration model '{}'", model);
    }
}

// ------------Molar Thermodynamic Properties -------------------------

CanteraDouble IdealSolnGasVPSS::enthalpy_mole() const
{
    updateStandardStateThermo();
    return RT() * mean_X(m_hss_RT);
}

CanteraDouble IdealSolnGasVPSS::entropy_mole() const
{
    updateStandardStateThermo();
    return GasConstant * (mean_X(m_sss_R) - sum_xlogx());
}

CanteraDouble IdealSolnGasVPSS::cp_mole() const
{
    updateStandardStateThermo();
    return GasConstant * mean_X(m_cpss_R);
}

CanteraDouble IdealSolnGasVPSS::cv_mole() const
{
    updateStandardStateThermo();
    CanteraDouble sum_xk_dVkdT = 0.0;
    CanteraDouble sum_xk_dVkdP = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble Xk = moleFraction(k);
        sum_xk_dVkdT += Xk * providePDSS(k)->dVdT();
        sum_xk_dVkdP += Xk * providePDSS(k)->dVdP();
    }
    if (sum_xk_dVkdP == 0.0) {
        return cp_mole();
    }
    // cv = cp - T*V*beta^2/kappa_T, using mixture beta and kappa_T from PDSS
    return cp_mole() - temperature() * sum_xk_dVkdT * sum_xk_dVkdT / (-sum_xk_dVkdP);
}

CanteraDouble IdealSolnGasVPSS::isothermalCompressibility() const
{
    updateStandardStateThermo();
    CanteraDouble sum_xk_dVkdP = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum_xk_dVkdP += moleFraction(k) * providePDSS(k)->dVdP();
    }
    return -sum_xk_dVkdP / molarVolume();
}

CanteraDouble IdealSolnGasVPSS::thermalExpansionCoeff() const
{
    updateStandardStateThermo();
    CanteraDouble sum_xk_dVkdT = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum_xk_dVkdT += moleFraction(k) * providePDSS(k)->dVdT();
    }
    return sum_xk_dVkdT / molarVolume();
}

void IdealSolnGasVPSS::setPressure(CanteraDouble p)
{
    m_Pcurrent = p;
    updateStandardStateThermo();
    calcDensity();
}

void IdealSolnGasVPSS::calcDensity()
{
    CanteraDouble v_mol = mean_X(getStandardVolumes());
    // Set the density in the parent object directly
    Phase::assignDensity(meanMolecularWeight() / v_mol);
}

Units IdealSolnGasVPSS::standardConcentrationUnits() const
{
    if (m_formGC != 0) {
        return Units(1.0, 0, -3, 0, 0, 0, 1);
    } else {
        return Units(1.0);
    }
}

void IdealSolnGasVPSS::getActivityConcentrations(span<CanteraDouble> c) const
{
    checkArraySize("IdealSolnGasVPSS::getActivityConcentrations", c.size(), m_kk);
    auto vss = getStandardVolumes();
    switch (m_formGC) {
    case 0:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = moleFraction(k);
        }
        break;
    case 1:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = moleFraction(k) / vss[k];
        }
        break;
    case 2:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = moleFraction(k) / vss[0];
        }
        break;
    }
}

CanteraDouble IdealSolnGasVPSS::standardConcentration(size_t k) const
{
    auto vss = getStandardVolumes();
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / vss[k];
    case 2:
        return 1.0/ vss[0];
    }
    return 0.0;
}

void IdealSolnGasVPSS::getActivityCoefficients(span<CanteraDouble> ac) const
{
    checkArraySize("IdealSolnGasVPSS::getActivityCoefficients", ac.size(), m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void IdealSolnGasVPSS::getChemPotentials(span<CanteraDouble> mu) const
{
    getStandardChemPotentials(mu);
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT() * log(xx);
    }
}

void IdealSolnGasVPSS::getPartialMolarEnthalpies(span<CanteraDouble> hbar) const
{
    getEnthalpy_RT(hbar);
    scale(hbar.begin(), hbar.end(), hbar.begin(), RT());
}

void IdealSolnGasVPSS::getPartialMolarEntropies(span<CanteraDouble> sbar) const
{
    getEntropy_R(sbar);
    scale(sbar.begin(), sbar.end(), sbar.begin(), GasConstant);
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (- log(xx));
    }
}

void IdealSolnGasVPSS::getPartialMolarIntEnergies(span<CanteraDouble> ubar) const
{
    getIntEnergy_RT(ubar);
    scale(ubar.begin(), ubar.end(), ubar.begin(), RT());
}

void IdealSolnGasVPSS::getPartialMolarCp(span<CanteraDouble> cpbar) const
{
    getCp_R(cpbar);
    scale(cpbar.begin(), cpbar.end(), cpbar.begin(), GasConstant);
}

void IdealSolnGasVPSS::getPartialMolarVolumes(span<CanteraDouble> vbar) const
{
    getStandardVolumes(vbar);
}

void IdealSolnGasVPSS::setToEquilState(span<const CanteraDouble> mu_RT)
{
    checkArraySize("IdealSolnGasVPSS::setToEquilState", mu_RT.size(), m_kk);
    updateStandardStateThermo();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    CanteraDouble pres = 0.0;
    CanteraDouble m_p0 = refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble tmp = -m_g0_RT[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            CanteraDouble tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setMoleFractions(m_pp);
    setPressure(pres);
}

bool IdealSolnGasVPSS::addSpecies(shared_ptr<Species> spec)
{
    bool added = VPStandardStateTP::addSpecies(spec);
    if (added) {
        m_pp.push_back(0.0);
    }
    return added;
}

void IdealSolnGasVPSS::initThermo()
{
    VPStandardStateTP::initThermo();
    if (m_input.hasKey("standard-concentration-basis")) {
        setStandardConcentrationModel(m_input["standard-concentration-basis"].asString());
    }
}

void IdealSolnGasVPSS::getParameters(AnyMap& phaseNode) const
{
    VPStandardStateTP::getParameters(phaseNode);
    // "unity" (m_formGC == 0) is the default, and can be omitted
    if (m_formGC == 1) {
        phaseNode["standard-concentration-basis"] = "species-molar-volume";
    } else if (m_formGC == 2) {
        phaseNode["standard-concentration-basis"] = "solvent-molar-volume";
    }
}

}
