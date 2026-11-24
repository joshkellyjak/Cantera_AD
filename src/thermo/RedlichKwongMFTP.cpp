//! @file RedlichKwongMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

#include <boost/algorithm/string.hpp>
#include <algorithm>

namespace Cantera
{

const CanteraDouble RedlichKwongMFTP::omega_a = 4.27480233540E-01;
const CanteraDouble RedlichKwongMFTP::omega_b = 8.66403499650E-02;
const CanteraDouble RedlichKwongMFTP::omega_vc = 3.33333333333333E-01;

RedlichKwongMFTP::RedlichKwongMFTP(const string& infile, const string& id_)
{
    initThermoFile(infile, id_);
}

void RedlichKwongMFTP::setSpeciesCoeffs(const string& species,
                                        CanteraDouble a0, CanteraDouble a1, CanteraDouble b)
{
    size_t k = speciesIndex(species, true);

    if (a1 != 0.0) {
        m_formTempParam = 1; // expression is temperature-dependent
    }

    m_a0(k, k) = a0;
    m_a1(k, k) = a1;

    // standard mixing rule for cross-species interaction term
    for (size_t j = 0; j < m_kk; j++) {
        if (k == j) {
            continue;
        }

        // m_a0 is initialized to NaN to mark uninitialized species
        if (isnan(m_a0(j, j))) {
            // The diagonal element of the jth species has not yet been defined.
            continue;
        } else if (isnan(m_a0(j, k))) {
            // Only use the mixing rules if the off-diagonal element has not already been defined by a
            // user-specified crossFluidParameters entry:
            m_a0(j, k) = m_a0(k, j) = sqrt(m_a0(j, j) * a0);
            m_a1(j, k) = m_a1(k, j) = sqrt(m_a1(j, j) * a1);
        }
    }
    m_a = m_a0;
    m_b[k] = b;
}

void RedlichKwongMFTP::setBinaryCoeffs(const string& species_i, const string& species_j,
                                       CanteraDouble a0, CanteraDouble a1)
{
    size_t ki = speciesIndex(species_i, true);
    size_t kj = speciesIndex(species_j, true);

    if (a1 != 0.0) {
        m_formTempParam = 1; // expression is temperature-dependent
    }

    m_binaryParameters[species_i][species_j] = {a0, a1};
    m_binaryParameters[species_j][species_i] = {a0, a1};
    m_a0(ki, kj) = m_a0(kj, ki) = a0;
    m_a1(ki, kj) = m_a1(kj, ki) = a1;
    m_a(ki, kj) = m_a(kj, ki) = a0;
}

// ------------Molar Thermodynamic Properties -------------------------

CanteraDouble RedlichKwongMFTP::cp_mole() const
{
    _updateReferenceStateThermo();
    CanteraDouble cpref = GasConstant * mean_X(m_cp0_R);
    if (m_bMix == 0.0) {
        return cpref;
    }
    CanteraDouble T = temperature();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble mv = molarVolume();
    CanteraDouble vpb = mv + m_bMix;
    pressureDerivatives();
    CanteraDouble dadt = da_dt();
    CanteraDouble fac = T * dadt - 3.0 * m_aMix / 2.0;
    CanteraDouble dHdT_V = cpref - 0.5 * dadt * log(vpb/mv) / (m_bMix * sqt)
        + mv * dpdT_ - GasConstant - log(vpb/mv) * fac / (2.0 * m_bMix * T * sqt);
    return dHdT_V - (mv + T * dpdT_ / dpdV_) * dpdT_;
}

CanteraDouble RedlichKwongMFTP::cv_mole() const
{
    _updateReferenceStateThermo();
    CanteraDouble cvref = GasConstant * (mean_X(m_cp0_R) - 1.0);
    if (m_bMix == 0.0) {
        return cvref;
    }
    CanteraDouble T = temperature();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble mv = molarVolume();
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble dadt = da_dt();
    CanteraDouble fac = T * dadt - 3.0 * m_aMix / 2.0;
    return (cvref - 1.0/(2.0 * m_bMix * T * sqt) * log(vpb/mv)*fac
            +1.0/(m_bMix * sqt) * log(vpb/mv)*(-0.5*dadt));
}

CanteraDouble RedlichKwongMFTP::pressure() const
{
    _updateReferenceStateThermo();

    //  Get a copy of the private variables stored in the State object
    CanteraDouble T = temperature();
    CanteraDouble molarV = meanMolecularWeight() / density();
    CanteraDouble pp = GasConstant * T/(molarV - m_bMix)
        - m_aMix/(sqrt(T) * molarV * (molarV + m_bMix));
    return pp;
}

CanteraDouble RedlichKwongMFTP::standardConcentration(size_t k) const
{
    getStandardVolumes(m_workS);
    return 1.0 / m_workS[k];
}

void RedlichKwongMFTP::getActivityCoefficients(span<CanteraDouble> ac) const
{
    checkArraySize("RedlichKwongMFTP::getActivityCoefficients", ac.size(), m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
    if (m_bMix == 0.0) {
        return;
    }
    CanteraDouble mv = molarVolume();
    CanteraDouble sqt = sqrt(temperature());
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble vmb = mv - m_bMix;
    CanteraDouble logv = log(vpb / mv);
    ArrayXd g = RT() * (log(RT() / (vmb * pressure())) + m_b / vmb)
        - 2.0 * m_Ak / (m_bMix * sqt) * logv
        + m_aMix * m_b / (m_bMix * m_bMix * sqt) * logv
        - m_aMix * m_b / (m_bMix * sqt * vpb);
    MappedVector(ac.data(), m_kk) = (g / RT()).exp();
}

// ---- Partial Molar Properties of the Solution -----------------

void RedlichKwongMFTP::getChemPotentials(span<CanteraDouble> mu) const
{
    getStandardChemPotentials(mu);
    Eigen::Map<ArrayXd> muVec(mu.data(), m_kk);
    auto x = asVectorXd(moleFractions_);
    muVec += RT() * x.array().max(SmallNumber).log();
    if (m_bMix == 0.0) {
        return;
    }

    CanteraDouble mv = molarVolume();
    CanteraDouble sqt = sqrt(temperature());
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble vmb = mv - m_bMix;
    CanteraDouble logv = log(vpb / mv);
    muVec += RT() * log(RT() / (vmb * pressure()))
        + RT() * m_b / vmb
        - 2.0 * m_Ak / (m_bMix * sqt) * logv
        + m_aMix * m_b / (m_bMix * m_bMix * sqt) * logv
        - m_aMix * m_b / (m_bMix * sqt * vpb);
}

void RedlichKwongMFTP::getPartialMolarEnthalpies(span<CanteraDouble> hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar.begin(), hbar.end(), hbar.begin(), RT());
    if (m_bMix == 0.0) {
        return;
    }

    // We calculate m_dpdni
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble vmb = mv - m_bMix;
    m_dpdni = RT() / vmb * (1.0 + m_b / vmb)
        - 2.0 * m_Ak / (sqt * mv * vpb)
        + m_aMix * m_b / (sqt * mv * vpb * vpb);
    CanteraDouble dadt = da_dt();
    CanteraDouble fac = T * dadt - 3.0 * m_aMix / 2.0;
    ArrayXd Sk = 2.0 * T * m_dAkdT - 3.0 * m_Ak;

    pressureDerivatives();
    CanteraDouble fac2 = mv + T * dpdT_ / dpdV_;
    CanteraDouble logv = log(vpb / mv);
    ArrayXd hE_v = mv * m_dpdni
        - RT()
        - m_b / (m_bMix * m_bMix * sqt) * logv * fac
        + (1.0 / (m_bMix * sqt) * logv) * Sk
        + m_b / (vpb * m_bMix * sqt) * fac;
    Eigen::Map<ArrayXd>(hbar.data(), m_kk) += hE_v - fac2 * m_dpdni;
}

void RedlichKwongMFTP::getPartialMolarEntropies(span<CanteraDouble> sbar) const
{
    getEntropy_R_ref(sbar);
    CanteraDouble T = temperature();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble mv = molarVolume();
    CanteraDouble pres = pressure();

    MappedVector sbarVec(sbar.data(), m_kk);
    sbarVec -= asVectorXd(moleFractions_).array().max(SmallNumber).log().matrix();
    sbarVec.array() -= log(pres / refPressure());
    sbarVec *= GasConstant;
    if (m_bMix == 0.0) {
        return;
    }
    CanteraDouble fac = da_dt() - m_aMix / (2.0 * T);
    CanteraDouble vmb = mv - m_bMix;
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble logv = log(vpb / mv);
    ArrayXd sdep = GasConstant * (log(RT() / (pres * vmb)) + 1)
        + GasConstant * m_b / vmb
        + m_Ak / (m_bMix * T * sqt) * logv
        - 2.0 * m_dAkdT / (m_bMix * sqt) * logv
        + m_b / (m_bMix * m_bMix * sqt) * logv * fac
        - m_b / (m_bMix * sqt * vpb) * fac;
    sbarVec -= sdep.matrix();

    pressureDerivatives();
    getPartialMolarVolumes(m_partialMolarVolumes);
    sbarVec += dpdT_ * asVectorXd(m_partialMolarVolumes);
}

void RedlichKwongMFTP::getPartialMolarIntEnergies(span<CanteraDouble> ubar) const
{
    // u_k = h_k - P * v_k
    getPartialMolarVolumes(m_partialMolarVolumes);
    getPartialMolarEnthalpies(ubar);
    CanteraDouble p = pressure();
    for (size_t k = 0; k < nSpecies(); k++) {
        ubar[k] -= p * m_partialMolarVolumes[k];
    }
}

void RedlichKwongMFTP::getPartialMolarIntEnergies_TV(span<CanteraDouble> utilde) const
{
    checkArraySize("RedlichKwongMFTP::getPartialMolarIntEnergies_TV", utilde.size(), m_kk);
    _updateReferenceStateThermo();

    Eigen::Map<ArrayXd> utildeVec(utilde.data(), m_kk);
    utildeVec = RT() * (asVectorXd(m_h0_RT).array() - 1.0);

    if (fabs(m_bMix) < SmallNumber) {
        return;
    }

    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble logv = log(vpb / mv);
    CanteraDouble F = T * da_dt() - 1.5 * m_aMix;
    CanteraDouble C = -logv / m_bMix + 1.0 / vpb;
    CanteraDouble pref = 1.0 / (m_bMix * sqrt(T));

    ArrayXd Sk = 2.0 * T * m_dAkdT - 3.0 * m_Ak;
    utildeVec += pref * (logv * Sk + m_b * F * C);
}

void RedlichKwongMFTP::getPartialMolarCp(span<CanteraDouble> cpbar) const
{
    checkArraySize("RedlichKwongMFTP::getPartialMolarCp", cpbar.size(), m_kk);
    _updateReferenceStateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] = GasConstant * m_cp0_R[k];
    }

    if (fabs(m_bMix) < SmallNumber) {
        return;
    }

    CanteraDouble T = temperature();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble v = molarVolume();
    CanteraDouble dadt = da_dt();
    CanteraDouble d2adt2 = 0.0; // current RK model uses a(T)=a0+a1*T

    // Mixture pressure derivatives
    pressureDerivatives();
    CanteraDouble dvdT_P = -dpdT_ / dpdV_;
    CanteraDouble vpb = v + m_bMix;
    CanteraDouble vmb = v - m_bMix;

    // P_TT, P_TV and P_VV for v''(T)|P from EOS:
    // v'' = -(P_TT + 2 P_TV v' + P_VV v'^2) / P_V
    CanteraDouble pTT = (-d2adt2 + dadt / T - 3.0 * m_aMix / (4.0 * T * T)) / (sqt * v * vpb);
    CanteraDouble pTV = -GasConstant / (vmb * vmb)
        + (dadt - m_aMix / (2.0 * T)) * (2.0 * v + m_bMix) / (sqt * v * v * vpb * vpb);
    CanteraDouble pVV = 2.0 * GasConstant * T / (vmb * vmb * vmb)
        - 2.0 * m_aMix * (3.0 * v * v + 3.0 * m_bMix * v + m_bMix * m_bMix)
            / (sqt * v * v * v * vpb * vpb * vpb);
    CanteraDouble d2vdT2_P = -(pTT + 2.0 * pTV * dvdT_P + pVV * dvdT_P * dvdT_P) / dpdV_;

    CanteraDouble F = T * dadt - 1.5 * m_aMix;
    CanteraDouble dF_dT = -0.5 * dadt + T * d2adt2;
    CanteraDouble logvpv = log(vpb / v);
    CanteraDouble termLPrime = -m_bMix * dvdT_P / (v * vpb);

    for (size_t k = 0; k < m_kk; k++) {
        CanteraDouble bk = m_b[k];
        CanteraDouble Ak = m_Ak[k];
        CanteraDouble dAk = m_dAkdT[k];
        CanteraDouble Sk = 2.0 * T * dAk - 3.0 * Ak;
        CanteraDouble dSk_dT = -dAk; // for linear a(T)

        // dp/dn_k at constant T,V,n_j
        CanteraDouble dpdni = RT() / vmb + RT() * bk / (vmb * vmb)
            - 2.0 * Ak / (sqt * v * vpb)
            + m_aMix * bk / (sqt * v * vpb * vpb);

        // d/dT|P [dp/dn_k]
        CanteraDouble d_dpdni_dT = GasConstant / vmb
            - GasConstant * T * dvdT_P / (vmb * vmb)
            + GasConstant * bk / (vmb * vmb)
            - 2.0 * GasConstant * T * bk * dvdT_P / (vmb * vmb * vmb)
            - 2.0 * dAk / (sqt * v * vpb)
            + Ak / (T * sqt * v * vpb)
            + 2.0 * Ak * (2.0 * v + m_bMix) * dvdT_P / (sqt * v * v * vpb * vpb)
            + bk * da_dt() / (sqt * v * vpb * vpb)
            - bk * m_aMix / (2.0 * T * sqt * v * vpb * vpb)
            - bk * m_aMix * (3.0 * v + m_bMix) * dvdT_P / (sqt * v * v * vpb * vpb * vpb);

        // d/dT|P of departure terms used in getPartialMolarEnthalpies
        CanteraDouble dterm1 = -bk / (m_bMix * m_bMix * sqt)
            * (termLPrime * F + logvpv * dF_dT - logvpv * F / (2.0 * T));
        CanteraDouble dterm2 = 1.0 / (m_bMix * sqt)
            * (termLPrime * Sk + logvpv * dSk_dT - logvpv * Sk / (2.0 * T));
        CanteraDouble dterm3 = bk / (m_bMix * sqt)
            * (dF_dT / vpb - F * dvdT_P / (vpb * vpb) - F / (2.0 * T * vpb));

        // cp_k = d hbar_k / dT at constant P and composition
        cpbar[k] += -GasConstant
            + (dvdT_P + T * d2vdT2_P) * dpdni
            + T * dvdT_P * d_dpdni_dT
            + dterm1 + dterm2 + dterm3;
    }
}

void RedlichKwongMFTP::getPartialMolarCv_TV(span<CanteraDouble> cvtilde) const
{
    checkArraySize("RedlichKwongMFTP::getPartialMolarCv_TV", cvtilde.size(), m_kk);
    _updateReferenceStateThermo();

    Eigen::Map<ArrayXd> cvtildeVec(cvtilde.data(), m_kk);
    cvtildeVec = GasConstant * (asVectorXd(m_cp0_R).array() - 1.0);

    if (fabs(m_bMix) < SmallNumber) {
        return;
    }

    CanteraDouble T = temperature();
    CanteraDouble sqt = sqrt(T);
    CanteraDouble mv = molarVolume();
    CanteraDouble dadt = da_dt();

    // Residual contribution for
    // \tilde{u}_k = (\partial U / \partial n_k)_{T,V,n_j}
    // with linear a(T): a_ij = a_ij,0 + a_ij,1 T
    CanteraDouble vpb = mv + m_bMix;
    CanteraDouble logv = log(vpb / mv);
    CanteraDouble F = T * dadt - 1.5 * m_aMix;
    CanteraDouble C = -logv / m_bMix + 1.0 / vpb;
    CanteraDouble pref = 1.0 / (m_bMix * sqt);

    ArrayXd Sk = 2.0 * T * m_dAkdT - 3.0 * m_Ak;
    ArrayXd ures = pref * (logv * Sk + m_b * F * C);
    ArrayXd dresdT = - pref * (logv * m_dAkdT + 0.5 * m_b * dadt * C)
        - 0.5 * ures / T;
    cvtildeVec += dresdT;
}

void RedlichKwongMFTP::getPartialMolarVolumes(span<CanteraDouble> vbar) const
{
    checkArraySize("RedlichKwongMFTP::getPartialMolarVolumes", vbar.size(), m_kk);
    CanteraDouble sqt = sqrt(temperature());
    CanteraDouble mv = molarVolume();
    CanteraDouble vmb = mv - m_bMix;
    CanteraDouble vpb = mv + m_bMix;
    ArrayXd num = RT() * (1.0 + (m_bMix + m_b) / vmb + m_bMix * m_b / (vmb * vmb))
        - 2.0 * m_Ak / (sqt * vpb)
        + m_aMix * m_b / (sqt * vpb * vpb);
    CanteraDouble denom = pressure() + RT() * m_bMix / (vmb * vmb) - m_aMix / (sqt * vpb * vpb);
    MappedVector(vbar.data(), m_kk) = num / denom;
}

bool RedlichKwongMFTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        size_t k = m_kk - 1;

        m_a.conservativeResize(m_kk, m_kk);
        m_a.row(k).setZero();
        m_a.col(k).setZero();

        // Initialize new entries to NaN to screen for species with undefined
        // pure-fluid parameters in the input file.
        m_b.conservativeResize(m_kk);
        m_b[k] = NAN;
        m_a0.conservativeResize(m_kk, m_kk);
        m_a0.row(k).setConstant(NAN);
        m_a0.col(k).setConstant(NAN);
        m_a1.conservativeResize(m_kk, m_kk);
        m_a1.row(k).setConstant(NAN);
        m_a1.col(k).setConstant(NAN);

        m_Ak.conservativeResizeLike(VectorXd::Zero(m_kk));
        m_dAkdT.conservativeResizeLike(VectorXd::Zero(m_kk));
        m_coeffSource.push_back(CoeffSource::EoS);
        m_partialMolarVolumes.push_back(0.0);
        m_dpdni.conservativeResizeLike(VectorXd::Zero(m_kk));
    }
    return added;
}

void RedlichKwongMFTP::initThermo()
{
    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<string, AnyMap*> dbSpecies;

    for (auto& [name, species] : m_species) {
        auto& data = species->input;
        size_t k = speciesIndex(name, true);
        if (!isnan(m_a0(k, k)) && !isnan(m_b[k])) {
            continue;
        }
        bool foundCoeffs = false;

        if (data.hasKey("equation-of-state") &&
            data["equation-of-state"].hasMapWhere("model", "Redlich-Kwong"))
        {
            // Read a and b coefficients from species 'input' information (that is, as
            // specified in a YAML input file) specific to the Redlich-Kwong model
            auto eos = data["equation-of-state"].getMapWhere(
                "model", "Redlich-Kwong");

            if (eos.hasKey("a") && eos.hasKey("b")) {
                CanteraDouble a0 = 0, a1 = 0;
                if (eos["a"].isScalar()) {
                    a0 = eos.convert("a", "Pa*m^6/kmol^2*K^0.5");
                } else {
                    auto avec = eos["a"].asVector<AnyValue>(2);
                    a0 = eos.units().convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                    a1 = eos.units().convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
                }
                CanteraDouble b = eos.convert("b", "m^3/kmol");
                foundCoeffs = true;
                setSpeciesCoeffs(name, a0, a1, b);
                m_coeffSource[k] = CoeffSource::EoS;
            }

            if (eos.hasKey("binary-a")) {
                AnyMap& binary_a = eos["binary-a"].as<AnyMap>();
                const UnitSystem& units = binary_a.units();
                for (auto& [name2, coeff] : binary_a) {
                    CanteraDouble a0 = 0, a1 = 0;
                    if (coeff.isScalar()) {
                        a0 = units.convert(coeff, "Pa*m^6/kmol^2*K^0.5");
                    } else {
                        auto avec = coeff.asVector<AnyValue>(2);
                        a0 = units.convert(avec[0], "Pa*m^6/kmol^2*K^0.5");
                        a1 = units.convert(avec[1], "Pa*m^6/kmol^2/K^0.5");
                    }
                    setBinaryCoeffs(name, name2, a0, a1);
                }
            }

            if (foundCoeffs) {
                continue;
            }
        }

        // Coefficients have not been populated from model-specific input
        CanteraDouble Tc = NAN, Pc = NAN;
        if (data.hasKey("critical-parameters")) {
            // Use critical state information stored in the species entry to
            // calculate a and b
            auto& critProps = data["critical-parameters"].as<AnyMap>();
            Tc = critProps.convert("critical-temperature", "K");
            Pc = critProps.convert("critical-pressure", "Pa");
            m_coeffSource[k] = CoeffSource::CritProps;
        } else {
            // Search 'crit-properties.yaml' to find Tc and Pc. Load data if needed
            if (critPropsDb.empty()) {
                critPropsDb = AnyMap::fromYamlFile("critical-properties.yaml");
                dbSpecies = critPropsDb["species"].asMap("name");
            }

            // All names in critical-properties.yaml are upper case
            auto ucName = boost::algorithm::to_upper_copy(name);
            if (dbSpecies.count(ucName)) {
                auto& spec = *dbSpecies.at(ucName);
                auto& critProps = spec["critical-parameters"].as<AnyMap>();
                Tc = critProps.convert("critical-temperature", "K");
                Pc = critProps.convert("critical-pressure", "Pa");
                m_coeffSource[k] = CoeffSource::Database;
            }
        }

        // Check if critical properties were found in either location
        if (!isnan(Tc)) {
            // Assuming no temperature dependence (that is, a1 = 0)
            CanteraDouble a = omega_a * pow(GasConstant, 2) * pow(Tc, 2.5) / Pc;
            CanteraDouble b = omega_b * GasConstant * Tc / Pc;
            setSpeciesCoeffs(name, a, 0.0, b);
        } else {
            throw InputFileError("RedlichKwongMFTP::initThermo", data,
                    "No critical property or Redlich-Kwong parameters found "
                    "for species {}.", name);
        }
    }
}

void RedlichKwongMFTP::getSpeciesParameters(const string& name,
                                            AnyMap& speciesNode) const
{
    MixtureFugacityTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name, true);
    if (m_coeffSource[k] == CoeffSource::EoS) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Redlich-Kwong", true);

        if (m_a1(k, k) != 0.0) {
            vector<AnyValue> coeffs(2);
            coeffs[0].setQuantity(m_a0(k, k), "Pa*m^6/kmol^2*K^0.5");
            coeffs[1].setQuantity(m_a1(k, k), "Pa*m^6/kmol^2/K^0.5");
            eosNode["a"] = std::move(coeffs);
        } else {
            eosNode["a"].setQuantity(m_a0(k, k),
                                    "Pa*m^6/kmol^2*K^0.5");
        }
        eosNode["b"].setQuantity(m_b[k], "m^3/kmol");
    } else if (m_coeffSource[k] == CoeffSource::CritProps) {
        auto& critProps = speciesNode["critical-parameters"];
        CanteraDouble a = m_a0(k, k);
        CanteraDouble b = m_b[k];
        CanteraDouble Tc = pow(a * omega_b / (b * omega_a * GasConstant), 2.0/3.0);
        CanteraDouble Pc = omega_b * GasConstant * Tc / b;
        critProps["critical-temperature"].setQuantity(Tc, "K");
        critProps["critical-pressure"].setQuantity(Pc, "Pa");
    }

    if (m_binaryParameters.count(name)) {
        auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
            "model", "Redlich-Kwong", true);
        AnyMap bin_a;
        for (const auto& [name2, coeffs] : m_binaryParameters.at(name)) {
            if (coeffs.second == 0) {
                bin_a[name2].setQuantity(coeffs.first, "Pa*m^6/kmol^2*K^0.5");
            } else {
                vector<AnyValue> C(2);
                C[0].setQuantity(coeffs.first, "Pa*m^6/kmol^2*K^0.5");
                C[1].setQuantity(coeffs.second, "Pa*m^6/kmol^2/K^0.5");
                bin_a[name2] = std::move(C);
            }
        }
        eosNode["binary-a"] = std::move(bin_a);
    }
}

CanteraDouble RedlichKwongMFTP::sresid() const
{
    if (m_bMix == 0.0) {
        return 0.0;
    }
    // note this agrees with tpx
    CanteraDouble hh = m_bMix / molarVolume();
    CanteraDouble zz = z();
    CanteraDouble T = temperature();
    CanteraDouble fac = da_dt() - m_aMix / (2.0 * T);
    CanteraDouble sresid_mol_R = log(zz*(1.0 - hh))
        + log(1.0 + hh) * fac / (sqrt(T) * GasConstant * m_bMix);
    return GasConstant * sresid_mol_R;
}

CanteraDouble RedlichKwongMFTP::hresid() const
{
    if (m_bMix == 0.0) {
        return 0.0;
    }
    // note this agrees with tpx
    CanteraDouble hh = m_bMix / molarVolume();
    CanteraDouble zz = z();
    CanteraDouble T = temperature();
    CanteraDouble fac = T * da_dt() - 3.0 *m_aMix / (2.0);
    return GasConstant * T * (zz - 1.0) + fac * log(1.0 + hh) / (sqrt(T) * m_bMix);
}

CanteraDouble RedlichKwongMFTP::liquidVolEst(CanteraDouble T, CanteraDouble& presGuess) const
{
    CanteraDouble v = m_bMix * 1.1;
    CanteraDouble atmp;
    CanteraDouble btmp;
    calculateAB(T, atmp, btmp);
    CanteraDouble pres = std::max(psatEst(T), presGuess);
    CanteraDouble Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
        int nsol = solveCubic(T, pres, atmp, btmp, Vroot);
        if (nsol == 1 || nsol == 2) {
            CanteraDouble pc = critPressure();
            if (pres > pc) {
                foundLiq = true;
            }
            pres *= 1.04;
        } else {
            foundLiq = true;
        }
    }

    if (foundLiq) {
        v = Vroot[0];
        presGuess = pres;
    } else {
        v = -1.0;
    }
    return v;
}

CanteraDouble RedlichKwongMFTP::densityCalc(CanteraDouble T, CanteraDouble presPa, int phaseRequested,
                                     CanteraDouble rhoguess)
{
    // It's necessary to set the temperature so that m_a_current is set correctly.
    setTemperature(T);
    CanteraDouble tcrit = critTemperature();
    CanteraDouble mmw = meanMolecularWeight();
    if (rhoguess == -1.0) {
        if (phaseRequested >= FLUID_LIQUID_0) {
                    CanteraDouble lqvol = liquidVolEst(T, presPa);
                    rhoguess = mmw / lqvol;
                }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to the routine
            rhoguess = presPa * mmw / (GasConstant * T);
    }


    CanteraDouble volguess = mmw / rhoguess;
    NSolns_ = solveCubic(T, presPa, m_aMix, m_bMix, Vroot_);

    // Ensure root ordering used for branch selection only contains physical roots.
    CanteraDouble vmin = std::max(0.0, m_bMix * (1.0 + 1e-12));
    vector<CanteraDouble> physicalRoots;
    for (CanteraDouble root : Vroot_) {
        if (std::isfinite(root) && root > vmin) {
            physicalRoots.push_back(root);
        }
    }
    std::sort(physicalRoots.begin(), physicalRoots.end());
    if (physicalRoots.empty()) {
        return -1.0;
    } else if (physicalRoots.size() == 1) {
        // Preserve branch-request semantics for single-root states:
        // return -2 when the requested phase is inconsistent with the
        // single branch identified by solveCubic() sign convention.
        if ((phaseRequested == FLUID_GAS && NSolns_ < 0)
            || (phaseRequested >= FLUID_LIQUID_0 && NSolns_ > 0))
        {
            return -2.0;
        }
        // Otherwise, accept the physically admissible root directly.
        return mmw / physicalRoots[0];
    } else if (physicalRoots.size() == 2) {
        Vroot_[0] = physicalRoots[0];
        Vroot_[1] = 0.5 * (physicalRoots[0] + physicalRoots[1]);
        Vroot_[2] = physicalRoots[1];
    } else {
        Vroot_[0] = physicalRoots[0];
        Vroot_[1] = physicalRoots[1];
        Vroot_[2] = physicalRoots[2];
    }

    CanteraDouble molarVolLast = Vroot_[0];
    if (NSolns_ >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = Vroot_[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[2];
        } else {
            if (volguess > Vroot_[1]) {
                molarVolLast = Vroot_[2];
            } else {
                molarVolLast = Vroot_[0];
            }
        }
    } else if (NSolns_ == 1) {
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT || phaseRequested == FLUID_UNDEFINED) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else if (NSolns_ == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[0];
        } else if (T > tcrit) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else {
        molarVolLast = Vroot_[0];
        return -1.0;
    }
    return mmw / molarVolLast;
}

CanteraDouble RedlichKwongMFTP::dpdVCalc(CanteraDouble T, CanteraDouble molarVol, CanteraDouble& presCalc) const
{
    CanteraDouble sqt = sqrt(T);
    presCalc = GasConstant * T / (molarVol - m_bMix)
               - m_aMix / (sqt * molarVol * (molarVol + m_bMix));

    CanteraDouble vpb = molarVol + m_bMix;
    CanteraDouble vmb = molarVol - m_bMix;
    return - GasConstant * T / (vmb * vmb)
        + m_aMix * (2 * molarVol + m_bMix) / (sqt * molarVol * molarVol * vpb * vpb);
}

CanteraDouble RedlichKwongMFTP::isothermalCompressibility() const
{
    pressureDerivatives();
    return -1 / (molarVolume() * dpdV_);
}

CanteraDouble RedlichKwongMFTP::thermalExpansionCoeff() const
{
    pressureDerivatives();
    return -dpdT_ / (molarVolume() * dpdV_);
}

CanteraDouble RedlichKwongMFTP::internalPressure() const
{
    pressureDerivatives();
    return temperature() * dpdT_ - pressure();
}

CanteraDouble RedlichKwongMFTP::soundSpeed() const
{
    pressureDerivatives();
    return molarVolume() * sqrt(-cp_mole() / cv_mole() * dpdV_ / meanMolecularWeight());
}

void RedlichKwongMFTP::pressureDerivatives() const
{
    CanteraDouble T = temperature();
    CanteraDouble mv = molarVolume();
    CanteraDouble pres;

    dpdV_ = dpdVCalc(T, mv, pres);
    CanteraDouble fac = da_dt() - m_aMix/(2.0 * T);
    dpdT_ = (GasConstant / (mv - m_bMix) - fac / (sqrt(T) * mv * (mv + m_bMix)));
}

void RedlichKwongMFTP::updateMixingExpressions()
{
    CanteraDouble T = temperature();
    auto x = asVectorXd(moleFractions_);
    if (m_formTempParam == 1) {
        m_a = m_a0 + T * m_a1;
        m_dAkdT = m_a1 * x;
    }

    m_bMix = m_b.matrix().dot(x);
    m_Ak = m_a * x;
    m_aMix = x.dot(m_Ak.matrix());
    if (isnan(m_bMix)) {
        // One or more species do not have specified coefficients.
        fmt::memory_buffer b;
        for (size_t k = 0; k < m_kk; k++) {
            if (isnan(m_b[k])) {
                if (b.size() > 0) {
                    fmt_append(b, ", {}", speciesName(k));
                } else {
                    fmt_append(b, "{}", speciesName(k));
                }
            }
        }
        throw CanteraError("RedlichKwongMFTP::updateMixingExpressions",
            "Missing Redlich-Kwong coefficients for species: {}", to_string(b));
    }
}

void RedlichKwongMFTP::calculateAB(CanteraDouble T, CanteraDouble& aCalc, CanteraDouble& bCalc) const
{
    auto x = asVectorXd(moleFractions_);
    bCalc = m_b.matrix().dot(x);
    if (m_formTempParam == 1) {
        aCalc = x.dot((m_a0 + T * m_a1) * x);
    } else {
        aCalc = x.dot(m_a0 * x);
    }
}

CanteraDouble RedlichKwongMFTP::da_dt() const
{
    if (m_formTempParam == 0) {
        return 0.0;
    }
    auto x = asVectorXd(moleFractions_);
    return x.dot(m_a1 * x);
}

void RedlichKwongMFTP::calcCriticalConditions(CanteraDouble& pc, CanteraDouble& tc, CanteraDouble& vc) const
{
    auto x = asVectorXd(moleFractions_);
    CanteraDouble a0 = x.dot(m_a0 * x);
    CanteraDouble aT = x.dot(m_a1 * x);
    CanteraDouble a = m_aMix;
    if (m_formTempParam != 0) {
        a = a0;
    }
    if (m_bMix <= 0.0) {
        tc = 1000000.;
        pc = 1.0E13;
        vc = omega_vc * GasConstant * tc / pc;
        return;
    }
    if (a <= 0.0) {
        tc = 0.0;
        pc = 0.0;
        vc = 2.0 * m_bMix;
        return;
    }
    CanteraDouble tmp = a * omega_b / (m_bMix * omega_a * GasConstant);
    CanteraDouble pp = 2./3.;
    CanteraDouble sqrttc, f, dfdt, deltatc;

    if (m_formTempParam == 0) {
        tc = pow(tmp, pp);
    } else {
        tc = pow(tmp, pp);
        for (int j = 0; j < 10; j++) {
            sqrttc = sqrt(tc);
            f = omega_a * m_bMix * GasConstant * tc * sqrttc / omega_b - aT * tc - a0;
            dfdt = 1.5 * omega_a * m_bMix * GasConstant * sqrttc / omega_b - aT;
            deltatc = - f / dfdt;
            tc += deltatc;
        }
        if (deltatc > 0.1) {
            throw CanteraError("RedlichKwongMFTP::calcCriticalConditions",
                "didn't converge");
        }
    }

    pc = omega_b * GasConstant * tc / m_bMix;
    vc = omega_vc * GasConstant * tc / pc;
}

int RedlichKwongMFTP::solveCubic(CanteraDouble T, CanteraDouble pres, CanteraDouble a, CanteraDouble b,
                                 span<CanteraDouble> Vroot) const
{

    // Derive the coefficients of the cubic polynomial to solve.
    CanteraDouble an = 1.0;
    CanteraDouble bn = - GasConstant * T / pres;
    CanteraDouble sqt = sqrt(T);
    CanteraDouble cn = - (GasConstant * T * b / pres - a/(pres * sqt) + b * b);
    CanteraDouble dn = - (a * b / (pres * sqt));

    CanteraDouble tmp = a * omega_b / (b * omega_a * GasConstant);
    CanteraDouble pp = 2./3.;
    CanteraDouble tc = pow(tmp, pp);
    CanteraDouble pc = omega_b * GasConstant * tc / b;
    CanteraDouble vc = omega_vc * GasConstant * tc / pc;

    int nSolnValues = MixtureFugacityTP::solveCubic(T, pres, a, b, a, Vroot, an, bn, cn, dn, tc, vc);

    return nSolnValues;
}

}
