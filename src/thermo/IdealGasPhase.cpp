/**
 *  @file IdealGasPhase.cpp
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see @ref thermoprops
 * and class @link Cantera::IdealGasPhase IdealGasPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

IdealGasPhase::IdealGasPhase(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);
}

// Molar Thermodynamic Properties of the Solution ------------------

CanteraDouble IdealGasPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx() - std::log(pressure() / refPressure()));
}

CanteraDouble IdealGasPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

CanteraDouble IdealGasPhase::cv_mole() const
{
    return cp_mole() - GasConstant;
}

CanteraDouble IdealGasPhase::soundSpeed() const {
    return sqrt(
        cp_mole() / cv_mole() * GasConstant / meanMolecularWeight() * temperature()
    );
}

CanteraDouble IdealGasPhase::standardConcentration(size_t k) const
{
    return pressure() / RT();
}

void IdealGasPhase::getActivityCoefficients(span<CanteraDouble> ac) const
{
    checkArraySize("IdealGasPhase::getActivityCoefficients", ac.size(), m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealGasPhase::getStandardChemPotentials(span<CanteraDouble> muStar) const
{
    getGibbs_ref(muStar);
    CanteraDouble tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

void IdealGasPhase::getChemPotentials(span<CanteraDouble> mu) const
{
    getStandardChemPotentials(mu);
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT() * log(xx);
    }
}

void IdealGasPhase::getPartialMolarEnthalpies(span<CanteraDouble> hbar) const
{
    checkArraySize("IdealGasPhase::getPartialMolarEnthalpies", hbar.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar.begin(), RT());
}

void IdealGasPhase::getPartialMolarEntropies(span<CanteraDouble> sbar) const
{
    checkArraySize("IdealGasPhase::getPartialMolarEntropies", sbar.size(), m_kk);
    auto _s = entropy_R_ref();
    scale(_s.begin(), _s.end(), sbar.begin(), GasConstant);
    CanteraDouble logp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (-logp - log(xx));
    }
}

void IdealGasPhase::getPartialMolarIntEnergies(span<CanteraDouble> ubar) const
{
    checkArraySize("IdealGasPhase::getPartialMolarIntEnergies", ubar.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
}

void IdealGasPhase::getPartialMolarIntEnergies_TV(span<CanteraDouble> utilde) const
{
    getPartialMolarIntEnergies(utilde);
}

void IdealGasPhase::getPartialMolarCp(span<CanteraDouble> cpbar) const
{
    checkArraySize("IdealGasPhase::getPartialMolarCp", cpbar.size(), m_kk);
    auto _cp = cp_R_ref();
    scale(_cp.begin(), _cp.end(), cpbar.begin(), GasConstant);
}

void IdealGasPhase::getPartialMolarCv_TV(span<CanteraDouble> cvtilde) const
{
    checkArraySize("IdealGasPhase::getPartialMolarCv_TV", cvtilde.size(), m_kk);
    auto _cp = cp_R_ref();
    for (size_t k = 0; k < m_kk; k++) {
        cvtilde[k] = GasConstant * (_cp[k] - 1.0);
    }
}

void IdealGasPhase::getPartialMolarVolumes(span<CanteraDouble> vbar) const
{
    checkArraySize("IdealGasPhase::getPartialMolarVolumes", vbar.size(), m_kk);
    CanteraDouble vol = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

void IdealGasPhase::getEnthalpy_RT(span<CanteraDouble> hrt) const
{
    checkArraySize("IdealGasPhase::getEnthalpy_RT", hrt.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt.begin());
}

void IdealGasPhase::getEntropy_R(span<CanteraDouble> sr) const
{
    checkArraySize("IdealGasPhase::getEntropy_R", sr.size(), m_kk);
    auto _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr.begin());
    CanteraDouble tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void IdealGasPhase::getGibbs_RT(span<CanteraDouble> grt) const
{
    checkArraySize("IdealGasPhase::getGibbs_RT", grt.size(), m_kk);
    auto gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt.begin());
    CanteraDouble tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void IdealGasPhase::getIntEnergy_RT(span<CanteraDouble> urt) const
{
    getIntEnergy_RT_ref(urt);
}

void IdealGasPhase::getCp_R(span<CanteraDouble> cpr) const
{
    checkArraySize("IdealGasPhase::getCp_R", cpr.size(), m_kk);
    auto _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr.begin());
}

void IdealGasPhase::getStandardVolumes(span<CanteraDouble> vol) const
{
    checkArraySize("IdealGasPhase::getStandardVolumes", vol.size(), m_kk);
    CanteraDouble tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

void IdealGasPhase::getEnthalpy_RT_ref(span<CanteraDouble> hrt) const
{
    checkArraySize("IdealGasPhase::getEnthalpy_RT_ref", hrt.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt.begin());
}

void IdealGasPhase::getGibbs_RT_ref(span<CanteraDouble> grt) const
{
    checkArraySize("IdealGasPhase::getGibbs_RT_ref", grt.size(), m_kk);
    auto gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt.begin());
}

void IdealGasPhase::getGibbs_ref(span<CanteraDouble> g) const
{
    checkArraySize("IdealGasPhase::getGibbs_ref", g.size(), m_kk);
    auto gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g.begin(), RT());
}

void IdealGasPhase::getEntropy_R_ref(span<CanteraDouble> er) const
{
    checkArraySize("IdealGasPhase::getEntropy_R_ref", er.size(), m_kk);
    auto _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er.begin());
}

void IdealGasPhase::getIntEnergy_RT_ref(span<CanteraDouble> urt) const
{
    checkArraySize("IdealGasPhase::getIntEnergy_RT_ref", urt.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void IdealGasPhase::getCp_R_ref(span<CanteraDouble> cprt) const
{
    checkArraySize("IdealGasPhase::getCp_R_ref", cprt.size(), m_kk);
    auto _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cprt.begin());
}

void IdealGasPhase::getStandardVolumes_ref(span<CanteraDouble> vol) const
{
    checkArraySize("IdealGasPhase::getStandardVolumes_ref", vol.size(), m_kk);
    CanteraDouble tmp = RT() / m_p0;
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

bool IdealGasPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            m_p0 = refPressure();
        }
        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_expg0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
        m_pp.push_back(0.0);
    }
    return added;
}

void IdealGasPhase::setToEquilState(span<const CanteraDouble> mu_RT)
{
    auto grt = gibbs_RT_ref();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    CanteraDouble pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 300.0) {
            CanteraDouble tmp2 = tmp / 300.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(300.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setMoleFractions(m_pp);
    setPressure(pres);
}

void IdealGasPhase::updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    CanteraDouble tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tnow) {
        m_spthermo.update(tnow, m_cp0_R, m_h0_RT, m_s0_R);
        cached.state1 = tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}
}
