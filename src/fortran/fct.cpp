/**
 *   Cantera Fortran interface library. This library of functions is designed
 *   to encapsulate Cantera functionality and make it available for
 *   use in languages and applications other than C++. A set of
 *   library functions is provided that are declared "extern C". All
 *   Cantera objects are stored and referenced by integers - no
 *   pointers are passed to or from the calling application.
 *
 *   @warning  The Fortran API is an experimental part of %Cantera and
 *      may be changed or removed without notice.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// Cantera includes
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/NoExitLogger.h"
#include "cantera/base/stringUtils.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include "clib_utils.h"
#include "cantera_clib/clib_defs.h"

using namespace Cantera;

typedef Cabinet<ThermoPhase> ThermoCabinet;
typedef Cabinet<Kinetics> KineticsCabinet;
typedef Cabinet<Transport> TransportCabinet;

typedef integer status_t;

template<> ThermoCabinet* ThermoCabinet::s_storage; // defined in ct.cpp
template<> KineticsCabinet* KineticsCabinet::s_storage; // defined in ct.cpp
template<> TransportCabinet* TransportCabinet::s_storage; // defined in ct.cpp

namespace {

ThermoPhase* _fph(const integer* n)
{
    return ThermoCabinet::at(*n).get();
}

static Kinetics* _fkin(const integer* n)
{
    return KineticsCabinet::at(*n).get();
}

ThermoPhase* _fth(const integer* n)
{
    return ThermoCabinet::at(*n).get();
}

shared_ptr<ThermoPhase> _fthermo(const integer* n)
{
    return ThermoCabinet::at(*n);
}

Transport* _ftrans(const integer* n)
{
    return TransportCabinet::at(*n).get();
}

} // unnamed namespace

string f2string(const char* s, ftnlen_t n)
{
    ftnlen_t k;
    string ss = "";
    for (k = 0; k < n; k++) {
        if (s[k] == '\0') {
            break;
        }
        ss += s[k];
    }
    return ss;
}

/**
 * Exported functions.
 */
extern "C" {

    status_t cantera_error_(const char* proc, const char* msg,
                            ftnlen_t proclen, ftnlen_t msglen)
    {
        string sproc = f2string(proc, proclen);
        string smsg = f2string(msg, msglen);
        throw CanteraError(sproc, smsg);
        return -1;
    }

    void ctturnoffexitonerror_() {
      setLogger(make_unique<Cantera::NoExitLogger>());
    }

    //--------------- Phase ---------------------//

    status_t phase_getname_(const integer* n, char* nm,
                            ftnlen_t lennm)
    {
        try {
            string pnm = _fph(n)->name();
            ftnlen_t lout = std::min(lennm, pnm.size());
            std::copy(pnm.c_str(), pnm.c_str() + lout, nm);
            for (ftnlen_t nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    integer phase_nelements_(const integer* n)
    {
        try {
            return _fph(n)->nElements();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer phase_nspecies_(const integer* n)
    {
        try {
            return _fph(n)->nSpecies();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble phase_temperature_(const integer* n)
    {
        try {
            return _fph(n)->temperature();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_settemperature_(const integer* n, CanteraDouble* t)
    {
        try {
            _fph(n)->setTemperature(*t);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble phase_density_(const integer* n)
    {
        try {
            return _fph(n)->density();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_setdensity_(const integer* n, CanteraDouble* rho)
    {
        try {
            _fph(n)->setDensity(*rho);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble phase_molardensity_(const integer* n)
    {
        try {
            return _fph(n)->molarDensity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble phase_meanmolecularweight_(const integer* n)
    {
        try {
            return _fph(n)->meanMolecularWeight();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    integer phase_elementindex_(const integer* n, char* nm, ftnlen_t lennm)
    {
        try {
            string elnm = f2string(nm, lennm);
            return _fph(n)->elementIndex(elnm) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer phase_speciesindex_(const integer* n, char* nm, ftnlen_t lennm)
    {
        try {
            string spnm = f2string(nm, lennm);
            return _fph(n)->speciesIndex(spnm) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t phase_getmolefractions_(const integer* n, CanteraDouble* x)
    {
        try {
            auto* phase = _fph(n);
            phase->getMoleFractions(span<CanteraDouble>(x, phase->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble phase_molefraction_(const integer* n, integer* k)
    {
        try {
            return _fph(n)->moleFraction(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }

    }

    status_t phase_getmassfractions_(const integer* n, CanteraDouble* y)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->getMassFractions(span<CanteraDouble>(y, p->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble phase_massfraction_(const integer* n, integer* k)
    {
        try {
            return _fph(n)->massFraction(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t phase_setmolefractions_(const integer* n, CanteraDouble* x, const integer* norm)
    {
        try {
            ThermoPhase* p = _fph(n);
            if (*norm) {
                p->setMoleFractions(span<const CanteraDouble>(x, p->nSpecies()));
            } else {
                p->setMoleFractions_NoNorm(span<const CanteraDouble>(x, p->nSpecies()));
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmolefractionsbyname_(const integer* n, char* x, ftnlen_t lx)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->setMoleFractionsByName(f2string(x, lx));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmassfractions_(const integer* n, CanteraDouble* y, const integer* norm)
    {
        try {
            ThermoPhase* p = _fph(n);
            if (*norm) {
                p->setMassFractions(span<const CanteraDouble>(y, p->nSpecies()));
            } else {
                p->setMassFractions_NoNorm(span<const CanteraDouble>(y, p->nSpecies()));
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_setmassfractionsbyname_(const integer* n, char* y, ftnlen_t leny)
    {
        try {
            ThermoPhase* p = _fph(n);
            p->setMassFractionsByName(f2string(y, leny));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getatomicweights_(const integer* n, CanteraDouble* atw)
    {
        try {
            ThermoPhase* p = _fph(n);
            auto wt = p->atomicWeights();
            copy(wt.begin(), wt.end(), atw);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getmolecularweights_(const integer* n, CanteraDouble* mw)
    {
        try {
            ThermoPhase* p = _fph(n);
            auto wt = p->molecularWeights();
            copy(wt.begin(), wt.end(), mw);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getspeciesname_(const integer* n, integer* k, char* nm, ftnlen_t lennm)
    {
        try {
            string spnm = _fph(n)->speciesName(*k-1);
            ftnlen_t lout = std::min(lennm, spnm.size());
            std::copy(spnm.c_str(), spnm.c_str() + lout, nm);
            for (ftnlen_t nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t phase_getelementname_(const integer* n, integer* m, char* nm, ftnlen_t lennm)
    {
        try {
            string elnm = _fph(n)->elementName(*m-1);
            ftnlen_t lout = std::min(lennm, elnm.size());
            std::copy(elnm.c_str(), elnm.c_str() + lout, nm);
            for (ftnlen_t nn = lout; nn < lennm; nn++) {
                nm[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble phase_natoms_(const integer* n, integer* k, integer* m)
    {
        try {
            return _fph(n)->nAtoms(*k-1,*m-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    //-------------- Thermo --------------------//

    integer th_newfromfile_(char* filename, char* phasename, ftnlen_t lenf, ftnlen_t lenp)
    {
        try {
            auto th = newThermo(f2string(filename, lenf),
                                f2string(phasename, lenp));
            return ThermoCabinet::add(th);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer th_geteostype_(const integer* n, char* buf, ftnlen_t lenbuf)
    {
        try {
            return copyString(_fth(n)->type(), buf, lenbuf);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble th_enthalpy_mole_(const integer* n)
    {
        try {
            return _fth(n)->enthalpy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_intenergy_mole_(const integer* n)
    {
        try {
            return _fth(n)->intEnergy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_entropy_mole_(const integer* n)
    {
        try {
            return _fth(n)->entropy_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_gibbs_mole_(const integer* n)
    {
        try {
            return _fth(n)->gibbs_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_cp_mole_(const integer* n)
    {
        try {
            return _fth(n)->cp_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_cv_mole_(const integer* n)
    {
        try {
            return _fth(n)->cv_mole();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_pressure_(const integer* n)
    {
        try {
            return _fth(n)->pressure();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_enthalpy_mass_(const integer* n)
    {
        try {
            return _fth(n)->enthalpy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_intenergy_mass_(const integer* n)
    {
        try {
            return _fth(n)->intEnergy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_entropy_mass_(const integer* n)
    {
        try {
            return _fth(n)->entropy_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_gibbs_mass_(const integer* n)
    {
        try {
            return _fth(n)->gibbs_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_cp_mass_(const integer* n)
    {
        try {
            return _fth(n)->cp_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_cv_mass_(const integer* n)
    {
        try {
            return _fth(n)->cv_mass();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t th_chempotentials_(const integer* n, CanteraDouble* murt)
    {
        try {
            ThermoPhase* thrm = _fth(n);
            thrm->getChemPotentials(span<CanteraDouble>(murt, thrm->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_setpressure_(const integer* n, CanteraDouble* p)
    {
        try {
            _fth(n)->setPressure(*p);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_hp_(const integer* n, CanteraDouble* v1, CanteraDouble* v2)
    {
        try {
            _fth(n)->setState_HP(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_uv_(const integer* n, CanteraDouble* v1, CanteraDouble* v2)
    {
        try {
            _fth(n)->setState_UV(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_sv_(const integer* n, CanteraDouble* v1, CanteraDouble* v2)
    {
        try {
            _fth(n)->setState_SV(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_set_sp_(const integer* n, CanteraDouble* v1, CanteraDouble* v2)
    {
        try {
            _fth(n)->setState_SP(*v1, *v2);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_equil_(const integer* n, char* XY, char* solver, CanteraDouble* rtol, int* max_steps, int* max_iter, int* estimate_equil, int* log_level, ftnlen_t lenxy, ftnlen_t lensolver)
    {
        try {
            _fth(n)->equilibrate(f2string(XY,lenxy), f2string(solver,lensolver), *rtol, *max_steps, *max_iter, *estimate_equil, *log_level);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble th_refpressure_(const integer* n)
    {
        try {
            return _fth(n)->refPressure();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble th_mintemp_(const integer* n, integer* k)
    {
        try {
            return _fth(n)->minTemp(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble th_maxtemp_(const integer* n, integer* k)
    {
        try {
            return _fth(n)->maxTemp(*k-1);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t th_getenthalpies_rt_(const integer* n, CanteraDouble* h_rt)
    {
        try {
            ThermoPhase* thrm = _fth(n);
            thrm->getEnthalpy_RT(span<CanteraDouble>(h_rt, thrm->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getgibbs_rt_(const integer* n, CanteraDouble* g_rt) {
      try {
        ThermoPhase* thrm = _fth(n);
        thrm->getGibbs_RT(span<CanteraDouble>(g_rt, thrm->nSpecies()));
      } catch(...) {
        return handleAllExceptions(-1, ERR);
      }
      return 0;
    }

    status_t th_getentropies_r_(const integer* n, CanteraDouble* s_r)
    {
        try {
            ThermoPhase* thrm = _fth(n);
            thrm->getEntropy_R(span<CanteraDouble>(s_r, thrm->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getcp_r_(const integer* n, integer* lenm, CanteraDouble* cp_r)
    {
        try {
            ThermoPhase* thrm = _fth(n);
            thrm->getCp_R(span<CanteraDouble>(cp_r, thrm->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getpartialmolarintenergies_r_(const integer* n, CanteraDouble* ie)
    {
        try {
            ThermoPhase* thrm = _fth(n);
            thrm->getPartialMolarIntEnergies(span<CanteraDouble>(ie, thrm->nSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t th_getpartialmolarenthalpies_(const integer* n, CanteraDouble* hbar) {
      try {
        ThermoPhase* thrm = _fth(n);
        thrm->getPartialMolarEnthalpies(span<CanteraDouble>(hbar, thrm->nSpecies()));
      } catch(...) {
        return handleAllExceptions(-1, ERR);
      }
      return 0;
    }

    status_t th_getpartialmolarcp_(const integer* n, CanteraDouble* cpbar) {
      try {
        ThermoPhase* thrm = _fth(n);
        thrm->getPartialMolarCp(span<CanteraDouble>(cpbar, thrm->nSpecies()));
      } catch(...) {
        return handleAllExceptions(-1, ERR);
      }
      return 0;
    }

    //-------------- Kinetics ------------------//

    integer kin_newfromfile_(const char* filename, const char* phasename,
                             integer* reactingPhase, const integer* neighbor1,
                             const integer* neighbor2, const integer* neighbor3,
                             const integer* neighbor4, ftnlen_t nlen, ftnlen_t plen)
    {
        try {
            vector<shared_ptr<ThermoPhase>> phases;
            phases.push_back(_fthermo(reactingPhase));
            if (*neighbor1 >= 0) {
                phases.push_back(_fthermo(neighbor1));
                if (*neighbor2 >= 0) {
                    phases.push_back(_fthermo(neighbor2));
                    if (*neighbor3 >= 0) {
                        phases.push_back(_fthermo(neighbor3));
                        if (*neighbor4 >= 0) {
                            phases.push_back(_fthermo(neighbor4));
                        }
                    }
                }
            }
            auto kin = newKinetics(phases, f2string(filename, nlen));
            return KineticsCabinet::add(kin);
        } catch (...) {
            return handleAllExceptions(999, ERR);
        }
    }

    //-------------------------------------
    integer kin_gettype_(const integer* n, char* buf, ftnlen_t buflen)
    {
        try {
            return copyString(_fkin(n)->kineticsType(), buf, buflen);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_start_(const integer* n, integer* p)
    {
        try {
            return _fkin(n)->kineticsSpeciesIndex(0, *p)+1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_speciesindex_(const integer* n, const char* nm, const char* ph,
                              ftnlen_t lennm, ftnlen_t lenph)
    {
        try {
            return _fkin(n)->kineticsSpeciesIndex(f2string(nm, lennm)) + 1;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //---------------------------------------

    integer kin_ntotalspecies_(const integer* n)
    {
        try {
            return _fkin(n)->nTotalSpecies();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_nreactions_(const integer* n)
    {
        try {
            return _fkin(n)->nReactions();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_nphases_(const integer* n)
    {
        try {
            return _fkin(n)->nPhases();
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer kin_phaseindex_(const integer* n, const char* ph, ftnlen_t lenph)
    {
        try {
            return _fkin(n)->phaseIndex(f2string(ph, lenph));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble kin_reactantstoichcoeff_(const integer* n, integer* k, integer* i)
    {
        try {
            return _fkin(n)->reactantStoichCoeff(*k-1,*i-1);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble kin_productstoichcoeff_(const integer* n, integer* k, integer* i)
    {
        try {
            return _fkin(n)->productStoichCoeff(*k-1,*i-1);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t kin_getreactiontype_(const integer* n, integer* i, char* buf, ftnlen_t lenbuf)
    {
        try {
            string r = _fkin(n)->reaction(*i-1)->type();
            ftnlen_t lout = std::min(lenbuf, r.size());
            std::copy(r.c_str(), r.c_str() + lout, buf);
            for (ftnlen_t nn = lout; nn < lenbuf; nn++) {
                buf[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getfwdratesofprogress_(const integer* n, CanteraDouble* fwdROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getFwdRatesOfProgress(span<CanteraDouble>(fwdROP, k->nReactions()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getrevratesofprogress_(const integer* n, CanteraDouble* revROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getRevRatesOfProgress(span<CanteraDouble>(revROP, k->nReactions()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    integer kin_isreversible_(const integer* n, integer* i)
    {
        try {
            return (int)_fkin(n)->isReversible(*i);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t kin_getnetratesofprogress_(const integer* n, CanteraDouble* netROP)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getNetRatesOfProgress(span<CanteraDouble>(netROP, k->nReactions()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getcreationrates_(const integer* n, CanteraDouble* cdot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getCreationRates(span<CanteraDouble>(cdot, k->nTotalSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getdestructionrates_(const integer* n, CanteraDouble* ddot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getDestructionRates(span<CanteraDouble>(ddot, k->nTotalSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getnetproductionrates_(const integer* n, CanteraDouble* wdot)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getNetProductionRates(span<CanteraDouble>(wdot, k->nTotalSpecies()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    CanteraDouble kin_multiplier_(const integer* n, integer* i)
    {
        try {
            return _fkin(n)->multiplier(*i);
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t kin_getequilibriumconstants_(const integer* n, CanteraDouble* kc)
    {
        try {
            Kinetics* k = _fkin(n);
            k->getEquilibriumConstants(span<CanteraDouble>(kc, k->nReactions()));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_getreactionstring_(const integer* n, integer* i, char* buf, ftnlen_t lenbuf)
    {
        try {
            Kinetics* k = _fkin(n);
            string r = k->reaction(*i-1)->equation();
            ftnlen_t lout = std::min(lenbuf, r.size());
            std::copy(r.c_str(), r.c_str() + lout, buf);
            for (ftnlen_t nn = lout; nn < lenbuf; nn++) {
                buf[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_setmultiplier_(const integer* n, integer* i, CanteraDouble* v)
    {
        try {
            _fkin(n)->setMultiplier(*i-1,*v);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t kin_advancecoverages_(const integer* n, CanteraDouble* tstep)
    {
        try {
            Kinetics* k = _fkin(n);
            if (k->kineticsType() == "surface" || k->kineticsType() == "edge") {
                ((InterfaceKinetics*)k)->advanceCoverages(*tstep);
            } else {
                throw CanteraError("kin_advanceCoverages",
                                   "wrong kinetics manager type");
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    //------------------- Transport ---------------------------

    integer newtransport_(char* model, integer* ith,
                          integer* loglevel, ftnlen_t lenmodel)
    {
        try {
            string mstr = f2string(model, lenmodel);
            auto t = _fthermo(ith);
            auto tr = newTransport(t, mstr);
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    integer trans_newdefault_(integer* ith, integer* loglevel, ftnlen_t lenmodel)
    {
        try {
            auto t = _fthermo(ith);
            auto tr = newTransport(t, "default");
            return TransportCabinet::add(tr);
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    CanteraDouble trans_viscosity_(const integer* n)
    {
        try {
            return _ftrans(n)->viscosity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble trans_electricalconductivity_(const integer* n)
    {
        try {
            return _ftrans(n)->electricalConductivity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    CanteraDouble trans_thermalconductivity_(const integer* n)
    {
        try {
            return _ftrans(n)->thermalConductivity();
        } catch (...) {
            return handleAllExceptions(DERR, DERR);
        }
    }

    status_t trans_getthermaldiffcoeffs_(const integer* n, CanteraDouble* dt)
    {
        try {
            auto tr = _ftrans(n);
            tr->getThermalDiffCoeffs(span<CanteraDouble>(dt, tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffs_(const integer* n, CanteraDouble* d)
    {
        try {
            auto tr = _ftrans(n);
            tr->getMixDiffCoeffs(span<CanteraDouble>(d, tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffsmass_(const integer* n, CanteraDouble* d)
    {
        try {
            auto tr = _ftrans(n);
            tr->getMixDiffCoeffsMass(span<CanteraDouble>(d, tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmixdiffcoeffsmole_(const integer* n, CanteraDouble* d)
    {
        try {
            auto tr = _ftrans(n);
            tr->getMixDiffCoeffsMole(span<CanteraDouble>(d, tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getbindiffcoeffs_(const integer* n, integer* ld, CanteraDouble* d)
    {
        try {
            auto tr = _ftrans(n);
            tr->getBinaryDiffCoeffs(*ld, span<CanteraDouble>(d, *ld * tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t trans_getmultidiffcoeffs_(const integer* n, integer* ld, CanteraDouble* d)
    {
        try {
            auto tr = _ftrans(n);
            tr->getMultiDiffCoeffs(*ld, span<CanteraDouble>(d, *ld * tr->thermo().nSpecies()));
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    //-------------------- Functions ---------------------------

    status_t ctphase_report_(const integer* nth,
                             char* buf, integer* show_thermo, ftnlen_t buflen)
    {
        try {
            bool stherm = (*show_thermo != 0);
            string s = _fth(nth)->report(stherm);
            if (s.size() >= buflen) {
                return -(s.size() + 1);
            }
            copy(s.begin(), s.end(), buf);
            for (ftnlen_t nn = s.size(); nn < buflen; nn++) {
                buf[nn] = ' ';
            }
            return 0;
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
    }

    status_t ctgetcanteraerror_(char* buf, ftnlen_t buflen)
    {
        try {
            string e;
            e = Application::Instance()->lastErrorMessage();
            ftnlen_t n = std::min(e.size(), buflen - 1);
            copy(e.begin(), e.begin() + n, buf);
            for (ftnlen_t nn = n; nn < buflen; nn++) {
                buf[nn] = ' ';
            }
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }

    status_t ctaddcanteradirectory_(integer* buflen, char* buf)
    {
        try {
            addDataDirectory(string(buf));
        } catch (...) {
            return handleAllExceptions(-1, ERR);
        }
        return 0;
    }
}
