#ifndef OPENMM_GKNPFORCEIMPL_H_
#define OPENMM_GKNPFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM-GKNP                              *
 * -------------------------------------------------------------------------- */

#include "openmm/AmoebaGKNPForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace OpenMM {

class System;

/**
 * This is the internal implementation of GKNPForce.
 */

class AmoebaGKNPForceImpl : public ForceImpl {
public:
    AmoebaGKNPForceImpl(const AmoebaGKNPForce& owner);
    ~AmoebaGKNPForceImpl();
    void initialize(ContextImpl& context);
    const AmoebaGKNPForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context,  bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
private:
    const AmoebaGKNPForce& owner;
    OpenMM::Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_GKNPFORCEIMPL_H_*/
