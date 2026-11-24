//! @file Integrators.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/numerics/IdasIntegrator.h"

namespace Cantera
{

Integrator* newIntegrator(const string& itype)
{
    if (itype == "CVODE") {
#if AD_ENABLED
        // AD: SUNDIALS (sunrealtype = double) cannot integrate AD types
        throw NotImplementedError("newIntegrator", "The CVODE integrator uses SUNDIALS, "
            "which only works with plain doubles, and is not available in AD builds "
            "of Cantera (ad_mode=forward/reverse).");
#else
        return new CVodesIntegrator();
#endif
    } else if (itype == "IDA") {
#if AD_ENABLED
        // AD: SUNDIALS (sunrealtype = double) cannot integrate AD types
        throw NotImplementedError("newIntegrator", "The IDA integrator uses SUNDIALS, "
            "which only works with plain doubles, and is not available in AD builds "
            "of Cantera (ad_mode=forward/reverse).");
#else
        return new IdasIntegrator();
#endif
    } else {
        throw CanteraError("newIntegrator",
                           "unknown integrator: "+itype);
    }
    return 0;
}

}
