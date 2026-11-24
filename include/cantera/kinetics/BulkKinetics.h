/**
 * @file BulkKinetics.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BULKKINETICS_H
#define CT_BULKKINETICS_H

#include "Kinetics.h"
#include "ThirdBodyCalc.h"

namespace Cantera
{

//! Specialization of Kinetics for chemistry in a single bulk phase
//! @ingroup kineticsmgr
class BulkKinetics : public Kinetics
{
public:
    //! @name Constructors and General Information
    //! @{
    BulkKinetics();

    string kineticsType() const override {
        return "bulk";
    }

    //! @}

    //! @name Reaction Mechanism Setup Routines
    //! @{
    bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;
    void addThirdBody(shared_ptr<Reaction> r);
    void modifyReaction(size_t i, shared_ptr<Reaction> rNew) override;
    void resizeSpecies() override;
    void resizeReactions() override;
    void setMultiplier(size_t i, CanteraDouble f) override;
    void invalidateCache() override;
    //! @}

    //! @name Reaction rate constants, rates of progress, and thermodynamic properties
    //! @{
    void getFwdRateConstants(span<CanteraDouble> kfwd) override;
    void getEquilibriumConstants(span<CanteraDouble> kc) override;
    void getRevRateConstants(span<CanteraDouble> krev, bool doIrreversible=false) override;

    void getDeltaGibbs(span<CanteraDouble> deltaG) override;
    void getDeltaEnthalpy(span<CanteraDouble> deltaH) override;
    void getDeltaEntropy(span<CanteraDouble> deltaS) override;

    void getDeltaSSGibbs(span<CanteraDouble> deltaG) override;
    void getDeltaSSEnthalpy(span<CanteraDouble> deltaH) override;
    void getDeltaSSEntropy(span<CanteraDouble> deltaS) override;
    //! @}

    //! @name Derivatives of rate constants and rates of progress
    //! @{
    void getDerivativeSettings(AnyMap& settings) const override;
    void setDerivativeSettings(const AnyMap& settings) override;
    void getFwdRateConstants_ddT(span<CanteraDouble> dkfwd) override;
    void getFwdRatesOfProgress_ddT(span<CanteraDouble> drop) override;
    void getRevRatesOfProgress_ddT(span<CanteraDouble> drop) override;
    void getNetRatesOfProgress_ddT(span<CanteraDouble> drop) override;
    void getFwdRateConstants_ddP(span<CanteraDouble> dkfwd) override;
    void getFwdRatesOfProgress_ddP(span<CanteraDouble> drop) override;
    void getRevRatesOfProgress_ddP(span<CanteraDouble> drop) override;
    void getNetRatesOfProgress_ddP(span<CanteraDouble> drop) override;
    void getFwdRateConstants_ddC(span<CanteraDouble> dkfwd) override;
    void getFwdRatesOfProgress_ddC(span<CanteraDouble> drop) override;
    void getRevRatesOfProgress_ddC(span<CanteraDouble> drop) override;
    void getNetRatesOfProgress_ddC(span<CanteraDouble> drop) override;
    Eigen::SparseMatrix<CanteraDouble> fwdRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> revRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> netRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<CanteraDouble> fwdRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<CanteraDouble> revRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<CanteraDouble> netRatesOfProgress_ddCi() override;
    void netProductionRates_ddCi(Eigen::SparseMatrix<CanteraDouble>& jac) override;
    //! @}

    //! @name Rate calculation intermediate methods
    //! @{

    void updateROP() override;

    void getThirdBodyConcentrations(span<CanteraDouble> concm) override;
    span<const CanteraDouble> thirdBodyConcentrations() const override {
        return m_concm;
    }

    //! @}

protected:
    //! @name Internal service methods
    //!
    //! @note These methods are for internal use, and seek to avoid code duplication
    //! while evaluating terms used for rate constants, rates of progress, and
    //! their derivatives.
    //! @{

    //! Multiply rate with third-body collider concentrations
    void processThirdBodies(span<CanteraDouble> rop);

    //! Multiply rate with inverse equilibrium constant
    void applyEquilibriumConstants(span<CanteraDouble> rop);

    //! Multiply rate with scaled temperature derivatives of the inverse
    //! equilibrium constant
    /*!
     *  This (scaled) derivative is handled by a finite difference.
     */
    void applyEquilibriumConstants_ddT(span<CanteraDouble> drkcn);

    //! Process temperature derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddT(const span<const CanteraDouble> in, span<CanteraDouble> drop);

    //! Process pressure derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddP(const span<const CanteraDouble> in, span<CanteraDouble> drop);

    //! Process concentration (molar density) derivative
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    //! @param mass_action  boolean indicating whether law of mass action applies
    void process_ddC(StoichManagerN& stoich, span<const CanteraDouble> in,
                     span<CanteraDouble> drop, bool mass_action=true);

    //! Process derivatives
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param ddX true: w.r.t mole fractions false: w.r.t species concentrations
    //! @return a sparse matrix of derivative contributions for each reaction of
    //! dimensions nTotalReactions by nTotalSpecies
    Eigen::SparseMatrix<CanteraDouble> calculateCompositionDerivatives(
        StoichManagerN& stoich, span<const CanteraDouble> in, bool ddX=true);

    //! Accumulate `sign * nu_net * d(rates-of-progress)/dC` into the
    //! already-patterned, compressed column-major sparse matrix `out`, mirroring
    //! calculateCompositionDerivatives() but scattering the stoichiometry-matrix
    //! product directly into `out` using a dense length-K workspace per column.
    //! `out`'s pattern must already cover every (k, m) the product can touch.
    void accumulateCompositionDerivatives(StoichManagerN& stoich,
                                          span<const CanteraDouble> in,
                                          Eigen::SparseMatrix<CanteraDouble>& out,
                                          CanteraDouble sign);

    //! Helper function ensuring that all rate derivatives can be calculated
    //! @param name  method name used for error output
    //! @throw CanteraError if ideal gas assumption does not hold
    void assertDerivativesValid(const string& name);

    //! @}

    //! Difference between the global reactants order and the global products
    //! order. Of type "double" to account for the fact that we can have real-
    //! valued stoichiometries.
    vector<CanteraDouble> m_dn;

    ThirdBodyCalc m_multi_concm; //!< used with MultiRate evaluator

    //! Third body concentrations
    vector<CanteraDouble> m_concm;

    //! Activity concentrations, as calculated by ThermoPhase::getActivityConcentrations
    vector<CanteraDouble> m_act_conc;

    //! Physical concentrations, as calculated by ThermoPhase::getConcentrations
    vector<CanteraDouble> m_phys_conc;

    //! Derivative settings
    bool m_jac_skip_third_bodies;
    bool m_jac_skip_falloff;
    bool m_jac_skip_nonideal;
    CanteraDouble m_jac_rtol_delta;

    bool m_ROP_ok = false;

    //! Buffers for partial rop results with length nReactions()
    vector<CanteraDouble> m_rbuf0;
    vector<CanteraDouble> m_rbuf1;
    vector<CanteraDouble> m_rbuf2;

    //! Dense scatter column of length nTotalSpecies() used by
    //! accumulateCompositionDerivatives() to fold the stoichiometry-matrix product
    //! into the result column-by-column.
    vector<CanteraDouble> m_rbuf3;
    vector<CanteraDouble> m_kf0; //!< Forward rate constants without perturbation
    vector<CanteraDouble> m_sbuf0;
    vector<CanteraDouble> m_state;
    vector<CanteraDouble> m_grt; //!< Standard chemical potentials for each species
};

}

#endif
