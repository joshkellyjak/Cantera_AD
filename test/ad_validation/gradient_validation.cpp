//! @file gradient_validation.cpp
//!
//! Gradient validation of an AD (CoDiPack) build of Cantera.
//!
//! The functions differentiated here mirror, call for call, what SU2's
//! CFluidCantera (branch feature_CANTERA_PR_AD) asks of Cantera, for the two SU2
//! combustion test cases:
//!   - hydrogen: TestCases/combustion/h2o2_hydrogen_air/2d_slit_burner_h2.cfg
//!               (h2o2.yaml, phase ohmech, mixture-averaged transport)
//!   - methane:  TestCases/combustion/1step_methane_air/2d_slit_burner_fluid.cfg
//!               (1S_CH4_MP1.yaml, phase CH4_MP1_mix, unity-Lewis-number transport)
//!
//! Two functions are checked, with SU2's inputs (N-1 transported mass fractions; the
//! last species is 1 - sum):
//!   F_T (SetTDState_T): (T, Y_1..Y_{N-1}) -> density, h, cp, cv, viscosity,
//!       conductivity, mixture-averaged D_k, h_k/RT, cp_k/R, chemical source terms
//!       W_k*wdot_k and the heat release rate
//!   F_h (SetTDState_h): (h, Y_1..Y_{N-1}) -> T and density, through SU2's Newton
//!       inversion of h(T)
//!
//! Checks:
//!   1. AD Jacobian vs. 4th-order central finite differences (FD), with a fresh
//!      Cantera Solution for every AD evaluation ("clean" protocol).
//!   2. Forward vs. reverse AD Jacobians (when given the other mode's output file).
//!   3. Known issue: Cantera caches properties keyed on the *value* of T. Re-using a
//!      Solution object (as SU2 does: one object for every cell) can therefore give
//!      wrong derivatives, and in reverse mode re-recording after a tape reset reads
//!      stale tape indices. Run last and reported separately; it does not affect the
//!      exit code.
//!
//! The AD mode is taken from the installed cantera/base/config.h. Build and run
//! with run_validation.sh.

#include "cantera/core.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/transport/Transport.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#if !AD_ENABLED
#error "Compile against an AD build of Cantera (ad_mode=forward or ad_mode=reverse)"
#endif

using namespace Cantera;
using std::string;
using std::vector;

#if defined(AD_FORWARD)
static const char* MODE = "forward";
#else
static const char* MODE = "reverse";
#endif

// Tolerance for AD vs. FD (scaled relative error, see scaledError())
static const double TOL_FD = 1e-6;
// Tolerance for forward vs. reverse AD
static const double TOL_FWD_REV = 1e-12;

static const double P_THERMO = 101325.0;  // THERMODYNAMIC_PRESSURE in both cases
static const double T_REF = 298.15;       // SU2's STD_REF_TEMP

struct CaseDef {
    string name;
    string mechanism;
    string phase;
    string transport;
    std::map<string, double> inletY;  // SU2 MARKER_INLET_SPECIES; remainder = last species
    double Tin;
};

// The parts of SU2's CFluidCantera that talk to Cantera
struct Model {
    std::shared_ptr<Solution> sol;
    ThermoPhase* th;
    Kinetics* kin;
    Transport* tr;
    size_t N;
    vector<CanteraDouble> W;      // molar masses
    vector<CanteraDouble> hForm;  // enthalpy of formation per unit mass (constant)

    explicit Model(const CaseDef& c) {
        sol = newSolution(c.mechanism, c.phase, c.transport);
        th = sol->thermo().get();
        kin = sol->kinetics().get();
        tr = sol->transport().get();
        N = th->nSpecies();
        W.resize(N);
        th->getMolecularWeights(W);
        // CFluidCantera::SetEnthalpyFormation
        vector<CanteraDouble> hrt(N);
        th->setState_TP(T_REF, P_THERMO);
        th->getEnthalpy_RT_ref(hrt);
        hForm.resize(N);
        for (size_t k = 0; k < N; k++) {
            hForm[k] = GasConstant * T_REF * hrt[k] / W[k];
        }
    }

    // CFluidCantera::SetMassFractions
    void setY(const CanteraDouble* yTransported) {
        vector<CanteraDouble> Y(N);
        CanteraDouble sum = 0.0;
        for (size_t k = 0; k < N - 1; k++) {
            Y[k] = yTransported[k];
            sum += yTransported[k];
        }
        Y[N - 1] = 1.0 - sum;
        th->setMassFractions(Y);
    }
};

struct Outputs {
    vector<CanteraDouble> v;
    vector<string> group;
    void add(const string& g, const CanteraDouble& x) {
        v.push_back(x);
        group.push_back(g);
    }
};

// F_T: CFluidCantera::SetTDState_T + ComputeChemicalSourceTerm + ComputeHeatRelease
static Outputs evalT(Model& m, const vector<CanteraDouble>& x)
{
    m.setY(&x[1]);
    m.th->setState_TP(x[0], P_THERMO);
    Outputs o;
    o.add("density", m.th->density());
    o.add("enthalpy", m.th->enthalpy_mass());
    o.add("cp", m.th->cp_mass());
    o.add("cv", m.th->cv_mass());
    o.add("viscosity", m.tr->viscosity());
    o.add("conductivity", m.tr->thermalConductivity());
    vector<CanteraDouble> D(m.N), hrt(m.N), cpr(m.N), wdot(m.N);
    m.tr->getMixDiffCoeffsMass(D);
    m.th->getEnthalpy_RT_ref(hrt);
    m.th->getCp_R_ref(cpr);
    m.kin->getNetProductionRates(wdot);
    CanteraDouble hrr = 0.0;
    for (size_t k = 0; k < m.N; k++) {
        o.add("D_k (mix-avg)", D[k]);
        o.add("h_k/RT", hrt[k]);
        o.add("cp_k/R", cpr[k]);
        if (k < m.N - 1) {
            o.add("W_k*wdot_k", m.W[k] * wdot[k]);
        }
        hrr += -1.0 * m.hForm[k] * m.W[k] * wdot[k];
    }
    o.add("heat release", hrr);
    return o;
}

// F_h: CFluidCantera::SetTDState_h (Newton inversion of h(T), as in SU2)
static Outputs evalH(Model& m, const vector<CanteraDouble>& x)
{
    const CanteraDouble enthalpy = x[0];
    const double toll = 1e-5;
    const int counter_limit = 20;
    CanteraDouble temp_iter = 300.0;
    CanteraDouble delta_temp_iter = 1e10;
    int counter = 0;
    m.setY(&x[1]);
    while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
        m.th->setState_TP(temp_iter, P_THERMO);
        CanteraDouble h_iter = m.th->enthalpy_mass();
        CanteraDouble cp_iter = m.th->cp_mass();
        delta_temp_iter = (enthalpy - h_iter) / cp_iter;
        temp_iter += delta_temp_iter;
    }
    m.th->setState_TP(temp_iter, P_THERMO);
    Outputs o;
    o.add("T(h,Y)", temp_iter);
    o.add("density(h,Y)", m.th->density());
    return o;
}

using Fn = std::function<Outputs(Model&, const vector<CanteraDouble>&)>;
using Matrix = vector<vector<double>>;  // [output][input]

static vector<CanteraDouble> toActive(const vector<double>& x)
{
    return vector<CanteraDouble>(x.begin(), x.end());
}

// Finite-difference step for mass fractions (O(1) variables)
static const double H_Y = 1e-5;

// 4th-order finite differences, passive evaluation. Input 0 (T or h) uses central
// differences with a step relative to scale0. Mass fractions use a fixed step,
// central where Y >= 2*H_Y and one-sided (forward) below that, because
// setMassFractions clips negative values: Y must not be pushed below zero. For a
// species with Y = 0 (e.g. Ar in the hydrogen case) this gives the one-sided
// derivative, which is also what AD returns there.
static Matrix fdJacobian(const CaseDef& c, Fn f, const vector<double>& x0, double scale0)
{
    Model m(c);
    size_t nIn = x0.size();
    size_t nOut = f(m, toActive(x0)).v.size();
    Matrix J(nOut, vector<double>(nIn));
    const vector<int> sC = {2, 1, -1, -2}, sF = {0, 1, 2, 3, 4};
    const vector<double> wC = {-1.0, 8.0, -8.0, 1.0}, wF = {-25.0, 48.0, -36.0, 16.0, -3.0};
    for (size_t j = 0; j < nIn; j++) {
        double h = (j == 0) ? 1e-4 * scale0 : H_Y;
        bool oneSided = (j > 0 && x0[j] < 2 * H_Y);
        const vector<int>& s = oneSided ? sF : sC;
        const vector<double>& w = oneSided ? wF : wC;
        vector<double> acc(nOut, 0.0);
        for (size_t p = 0; p < s.size(); p++) {
            vector<double> x = x0;
            x[j] += s[p] * h;
            Outputs o = f(m, toActive(x));
            for (size_t i = 0; i < nOut; i++) {
                acc[i] += w[p] * o.v[i].getValue();
            }
        }
        for (size_t i = 0; i < nOut; i++) {
            J[i][j] = acc[i] / (12.0 * h);
        }
    }
    return J;
}

// AD Jacobian. With reuse == nullptr every evaluation gets a fresh Solution
// ("clean"); otherwise the given object is used for everything, as SU2 does.
static Matrix adJacobian(const CaseDef& c, Fn f, const vector<double>& x0, Model* reuse)
{
    size_t nIn = x0.size();
    Matrix J;
#if defined(AD_FORWARD)
    for (size_t j = 0; j < nIn; j++) {
        std::unique_ptr<Model> fresh;
        Model* m = reuse;
        if (!m) {
            fresh = std::make_unique<Model>(c);
            m = fresh.get();
        }
        vector<CanteraDouble> x = toActive(x0);
        x[j].setGradient(1.0);
        Outputs o = f(*m, x);
        if (J.empty()) {
            J.assign(o.v.size(), vector<double>(nIn));
        }
        for (size_t i = 0; i < o.v.size(); i++) {
            J[i][j] = o.v[i].getGradient();
        }
    }
#else
    auto& tape = CanteraDouble::getTape();
    std::unique_ptr<Model> fresh;
    Model* m = reuse;
    if (!m) {
        fresh = std::make_unique<Model>(c);
        m = fresh.get();
    }
    int recordings = reuse ? 2 : 1;  // reuse: record once, reset, record again
    for (int r = 0; r < recordings; r++) {
        tape.reset();
        tape.setActive();
        vector<CanteraDouble> x = toActive(x0);
        for (auto& xi : x) {
            tape.registerInput(xi);
        }
        Outputs o = f(*m, x);
        for (auto& yi : o.v) {
            tape.registerOutput(yi);
        }
        tape.setPassive();
        if (r < recordings - 1) {
            continue;
        }
        J.assign(o.v.size(), vector<double>(nIn));
        for (size_t i = 0; i < o.v.size(); i++) {
            tape.clearAdjoints();
            o.v[i].setGradient(1.0);
            tape.evaluate();
            for (size_t j = 0; j < nIn; j++) {
                J[i][j] = x[j].getGradient();
            }
        }
    }
    tape.reset();
#endif
    return J;
}

// Error of row A relative to reference row B (one output y_i, all inputs). Entries
// are multiplied by an input scale xs_j (T or c_p*T for the first input, 1 for mass
// fractions) so that inputs of different units are comparable. The error of each
// entry is relative to max(|B_ij|, 1e-3 * max_j |B_ij|, 1e-4 * |y_i|): sensitivities
// that are tiny compared with the row's largest one, or with the output itself, are
// measured against that floor instead (FD roundoff is ~1e-11 * |y_i| there).
static double scaledError(const vector<double>& a, const vector<double>& b,
                          const vector<double>& xs, double yAbs)
{
    double rowMax = 0.0;
    for (size_t j = 0; j < b.size(); j++) {
        rowMax = std::max(rowMax, std::abs(b[j] * xs[j]));
    }
    double err = 0.0;
    for (size_t j = 0; j < b.size(); j++) {
        double denom = std::max({std::abs(b[j] * xs[j]), 1e-3 * rowMax, 1e-4 * yAbs});
        double diff = std::abs(a[j] - b[j]) * xs[j];
        if (diff > 0) {
            err = std::max(err, denom > 0 ? diff / denom : HUGE_VAL);
        }
    }
    return err;
}

struct State {
    string label;
    double T;
    vector<double> Ytr;  // transported mass fractions (N-1)
};

// States between the unburnt inlet mixture and its adiabatic equilibrium (HP),
// so that every species is present
static vector<State> makeStates(const CaseDef& c)
{
    Model m(c);
    size_t N = m.N;
    vector<CanteraDouble> Yu(N, 0.0);
    double sum = 0.0;
    for (auto& [name, y] : c.inletY) {
        Yu[m.th->speciesIndex(name)] = y;
        sum += y;
    }
    Yu[N - 1] = 1.0 - sum;
    m.th->setMassFractions(Yu);
    m.th->setState_TP(c.Tin, P_THERMO);
    m.th->equilibrate("HP");
    double Tb = m.th->temperature().getValue();
    vector<CanteraDouble> Yb(N);
    m.th->getMassFractions(Yb);
    vector<State> states;
    for (double prog : {0.05, 0.5, 0.95}) {
        State s;
        std::ostringstream lab;
        lab << "progress " << prog;
        s.label = lab.str();
        s.T = (1 - prog) * c.Tin + prog * Tb;
        for (size_t k = 0; k < N - 1; k++) {
            s.Ytr.push_back((1 - prog) * Yu[k].getValue() + prog * Yb[k].getValue());
        }
        states.push_back(s);
    }
    return states;
}

static std::map<string, double> loadJacobians(const string& file)
{
    std::map<string, double> J;
    std::ifstream in(file);
    string key;
    double v;
    while (in >> key >> v) {
        J[key] = v;
    }
    return J;
}

static string jacKey(const string& cs, size_t st, const string& fn, size_t i, size_t j)
{
    std::ostringstream k;
    k << cs << '|' << st << '|' << fn << '|' << i << '|' << j;
    return k.str();
}

static vector<double> stateInputs(const State& s)
{
    vector<double> x = {s.T};
    x.insert(x.end(), s.Ytr.begin(), s.Ytr.end());
    return x;
}

// A species with Y = 0 exactly (Ar in the hydrogen case: neither in the inlet nor
// produced) sits on the kink of the clipping in Phase::setMassFractions
// (std::max(y, 0)), where the derivative is one-sided. The main check lifts such
// species to Y_LIFT; the kink itself is reported by zeroSpeciesCheck().
static const double Y_LIFT = 1e-6;

static State lifted(State s, int& nLifted)
{
    nLifted = 0;
    for (auto& y : s.Ytr) {
        if (y == 0.0) {
            y = Y_LIFT;
            nLifted++;
        }
    }
    return s;
}

// AD derivative w.r.t. a species at Y = 0 vs. the one-sided (Y >= 0) FD derivative
static void zeroSpeciesCheck(const CaseDef& c)
{
    vector<double> x = stateInputs(makeStates(c)[1]);
    Model probe(c);
    bool any = false;
    Matrix Jad, Jfd;
    Outputs yv;
    for (size_t j = 1; j < x.size(); j++) {
        if (x[j] != 0.0) {
            continue;
        }
        if (!any) {
            Jad = adJacobian(c, evalT, x, nullptr);
            Jfd = fdJacobian(c, evalT, x, x[0]);
            yv = evalT(probe, toActive(x));
            any = true;
        }
        vector<double> xs(x.size(), 1.0);
        xs[0] = x[0];
        double e = 0.0;
        for (size_t i = 0; i < Jad.size(); i++) {
            double rowMax = 0.0;
            for (size_t jj = 0; jj < x.size(); jj++) {
                rowMax = std::max(rowMax, std::abs(Jfd[i][jj] * xs[jj]));
            }
            double denom = std::max({std::abs(Jfd[i][j]), 1e-3 * rowMax,
                                     1e-4 * std::abs(yv.v[i].getValue())});
            if (denom > 0) {
                e = std::max(e, std::abs(Jad[i][j] - Jfd[i][j]) / denom);
            }
        }
        printf("  %-8s Y_%s = 0: AD vs one-sided FD, max err %.2e  %s\n", c.name.c_str(),
               probe.th->speciesName(j - 1).c_str(), e,
               e < TOL_FD ? "OK" : "KNOWN ISSUE (clipping kink)");
    }
    if (!any) {
        printf("  %-8s no transported species at Y = 0\n", c.name.c_str());
    }
}

// Two "cells" with the same temperature, evaluated one after the other with the
// same Solution, as SU2 does for neighbouring cells in a uniform region. The
// second evaluation hits the T-keyed caches filled by the first, so dh_b/dT_b
// misses the species-enthalpy part (and in reverse mode it goes to T_a instead).
static void twoCellCheck(const CaseDef& c)
{
    vector<double> xT = stateInputs(makeStates(c)[1]);
    Model m(c);
    vector<CanteraDouble> xa = toActive(xT), xb = toActive(xT);
    double got, leak = 0.0;
#if defined(AD_FORWARD)
    xb[0].setGradient(1.0);
    evalT(m, xa);
    got = evalT(m, xb).v[1].getGradient();
#else
    auto& tape = CanteraDouble::getTape();
    tape.reset();
    tape.setActive();
    tape.registerInput(xa[0]);
    tape.registerInput(xb[0]);
    evalT(m, xa);
    CanteraDouble hb = evalT(m, xb).v[1];
    tape.registerOutput(hb);
    tape.setPassive();
    hb.setGradient(1.0);
    tape.evaluate();
    got = xb[0].getGradient();
    leak = xa[0].getGradient();
    tape.reset();
#endif
    Model ref(c);
    double cp = evalT(ref, toActive(xT)).v[2].getValue();
    bool ok = std::abs(got - cp) / cp < TOL_FD && leak == 0.0;
    printf("  %-8s two cells, equal T: dh_b/dT_b = % .6e (expected cp = %.6e)", c.name.c_str(),
           got, cp);
#if !defined(AD_FORWARD)
    printf(", dh_b/dT_a = % .6e (expected 0)", leak);
#endif
    printf("  %s\n", ok ? "OK" : "KNOWN ISSUE");
}

// Full F_T Jacobian with one Solution re-used for every AD evaluation. Forward
// mode: the same object is re-seeded for each input. Reverse mode: the tape is
// recorded, reset and recorded again, so cached values keep indices from the
// discarded recording (undefined behaviour: stale indices can exceed the new tape).
static void reuseCheck(const CaseDef& c)
{
    vector<double> xT = stateInputs(makeStates(c)[1]);
    vector<double> xs(xT.size(), 1.0);
    xs[0] = xT[0];
    Matrix Jclean = adJacobian(c, evalT, xT, nullptr);
    Model mo(c);
    Outputs o = evalT(mo, toActive(xT));
    Model reused(c);
    Matrix Jr = adJacobian(c, evalT, xT, &reused);
    double e = 0.0;
    for (size_t i = 0; i < Jr.size(); i++) {
        e = std::max(e, scaledError(Jr[i], Jclean[i], xs, std::abs(o.v[i].getValue())));
    }
#if !defined(AD_FORWARD)
    CanteraDouble::getTape().deleteAdjointVector();
#endif
    printf("  %-8s F_T Jacobian, one Solution re-used (%s): max err vs clean AD %.2e  %s\n",
           c.name.c_str(), strcmp(MODE, "forward") ? "re-recorded after tape reset"
                                                     : "re-seeded",
           e, e < TOL_FD ? "OK" : "KNOWN ISSUE");
}

int main(int argc, char** argv)
{
    string writeFile, compareFile;
    for (int a = 1; a < argc - 1; a++) {
        if (!strcmp(argv[a], "--write")) writeFile = argv[a + 1];
        if (!strcmp(argv[a], "--compare")) compareFile = argv[a + 1];
    }
    std::ofstream jacOut;
    if (!writeFile.empty()) {
        jacOut.open(writeFile);
        jacOut.precision(17);
    }
    std::map<string, double> other;
    if (!compareFile.empty()) {
        other = loadJacobians(compareFile);
        printf("Loaded %zu Jacobian entries from %s\n", other.size(), compareFile.c_str());
    }

    printf("Cantera %s, AD mode: %s, CoDiPack %d.%d\n", CANTERA_VERSION, MODE,
           CODI_MAJOR_VERSION, CODI_MINOR_VERSION);
    printf("AD vs FD tolerance %.0e (scaled relative error)\n\n", TOL_FD);

    vector<CaseDef> cases = {
        {"hydrogen", "h2o2.yaml", "ohmech", "mixture-averaged",
         {{"H2", 0.017311}, {"O2", 0.228966}}, 300.0},
        {"methane", "1S_CH4_MP1.yaml", "CH4_MP1_mix", "unity-Lewis-number",
         {{"O2", 0.222598}, {"CH4", 0.044642}}, 300.0},
    };

    int failures = 0;
    double worstFwdRev = 0.0;
    size_t comparedFwdRev = 0;
    for (auto& c : cases) {
        Model probe(c);
        printf("=== %s: %s / %s / %s, %zu species (last = %s)\n", c.name.c_str(),
               c.mechanism.c_str(), c.phase.c_str(), c.transport.c_str(), probe.N,
               probe.th->speciesName(probe.N - 1).c_str());
        auto states = makeStates(c);
        for (size_t st = 0; st < states.size(); st++) {
            int nLifted;
            const State s = lifted(states[st], nLifted);
            // F_T inputs and F_h inputs
            vector<double> xT = stateInputs(s);
            Model mh(c);
            Outputs oT = evalT(mh, toActive(xT));
            double hval = oT.v[1].getValue();
            double cpT = oT.v[2].getValue() * s.T;
            vector<double> xH = {hval};
            xH.insert(xH.end(), s.Ytr.begin(), s.Ytr.end());
            double ymin = *std::min_element(s.Ytr.begin(), s.Ytr.end());
            printf("\n--- %s: T = %.1f K, h = %.4e J/kg, min transported Y = %.2e",
                   s.label.c_str(), s.T, hval, ymin);
            if (nLifted) {
                printf(" (%d species lifted from 0 to %.0e)", nLifted, Y_LIFT);
            }
            printf("\n");

            struct Job { string fn; Fn f; vector<double> x; double scale0; };
            vector<Job> jobs = {{"F_T", evalT, xT, s.T}, {"F_h", evalH, xH, cpT}};
            for (auto& job : jobs) {
                // input scales for the error measure: T (or c_p T for h), 1 for Y
                vector<double> xs(job.x.size(), 1.0);
                xs[0] = job.scale0;
                Matrix Jfd = fdJacobian(c, job.f, job.x, job.scale0);
                Matrix Jad = adJacobian(c, job.f, job.x, nullptr);
                Model dummy(c);
                Outputs o = job.f(dummy, toActive(job.x));
                // per-group worst error
                std::map<string, double> worst;
                vector<string> order;
                for (size_t i = 0; i < Jad.size(); i++) {
                    if (!worst.count(o.group[i])) order.push_back(o.group[i]);
                    double e = scaledError(Jad[i], Jfd[i], xs, std::abs(o.v[i].getValue()));
                    worst[o.group[i]] = std::max(worst[o.group[i]], e);
                    vector<double> refRow;
                    for (size_t j = 0; j < Jad[i].size(); j++) {
                        string key = jacKey(c.name, st, job.fn, i, j);
                        if (jacOut.is_open()) jacOut << key << ' ' << Jad[i][j] << '\n';
                        if (other.count(key)) refRow.push_back(other[key]);
                    }
                    if (refRow.size() == Jad[i].size()) {
                        worstFwdRev = std::max(worstFwdRev, scaledError(Jad[i], refRow, xs,
                                                                   std::abs(o.v[i].getValue())));
                        comparedFwdRev += refRow.size();
                    }
                }
                printf("  %s  AD vs FD (%zu outputs x %zu inputs)\n", job.fn.c_str(),
                       Jad.size(), job.x.size());
                for (auto& g : order) {
                    bool ok = worst[g] < TOL_FD;
                    failures += !ok;
                    printf("    %-15s max err %.2e  %s\n", g.c_str(), worst[g],
                           ok ? "OK" : "FAIL");
                }
                // a couple of raw numbers for eyeballing
                size_t last = Jad.size() - 1;
                printf("    e.g. d(%s)/d(x0): AD % .10e  FD % .10e\n", o.group[0].c_str(),
                       Jad[0][0], Jfd[0][0]);
                printf("    e.g. d(%s)/d(x0): AD % .10e  FD % .10e\n",
                       o.group[last].c_str(), Jad[last][0], Jfd[last][0]);
            }
        }
        printf("\n");
    }

    if (!compareFile.empty()) {
        bool ok = comparedFwdRev > 0 && worstFwdRev < TOL_FWD_REV;
        failures += !ok;
        printf("Forward vs reverse AD: %zu entries compared, max scaled difference %.2e  %s\n",
               comparedFwdRev, worstFwdRev, ok ? "OK" : "FAIL");
    }

    // Known issue, run last: in reverse mode the re-use check reads stale tape
    // indices, which can corrupt later recordings.
    printf("\n=== Known issue: value-keyed caches when one Solution is re-used "
           "(as SU2 does) ===\n");
    for (auto& c : cases) {
        twoCellCheck(c);
    }
    printf("\n=== Known issue: derivative w.r.t. a species at Y = 0 (setMassFractions "
           "clips with std::max(y, 0), slope 0 at y = 0) ===\n");
    for (auto& c : cases) {
        zeroSpeciesCheck(c);
    }
    printf("\n=== Known issue (run last): full Jacobian with one re-used Solution ===\n");
    for (auto& c : cases) {
        reuseCheck(c);
    }
    printf("\n");
    printf("%s: %d failure(s) in AD vs FD%s\n", failures ? "FAILED" : "PASSED", failures,
           compareFile.empty() ? "" : " / forward vs reverse");
    return failures ? 1 : 0;
}
