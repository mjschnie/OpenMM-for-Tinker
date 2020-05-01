/* -------------------------------------------------------------------------- *
 *                       OpenMM-GKCavitation                                  *
 * -------------------------------------------------------------------------- */

#include <iostream>
#include "openmm/AmoebaGKCavitationForce.h"
#include "openmm/internal/AmoebaGKCavitationForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace OpenMM;
using namespace std;

AmoebaGKCavitationForce::AmoebaGKCavitationForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0) {
}

int AmoebaGKCavitationForce::addParticle(double radius, double gamma, double vdw_alpha, double charge, bool ishydrogen){
  ParticleInfo particle(radius, gamma, vdw_alpha, charge, ishydrogen);
  particles.push_back(particle);
  return particles.size()-1;
}

void AmoebaGKCavitationForce::setParticleParameters(int index, double radius, double gamma, double vdw_alpha, double charge, bool ishydrogen){
  particles[index].radius = radius;
  particles[index].gamma = gamma;
  particles[index].vdw_alpha = vdw_alpha;
  particles[index].charge = charge;
  particles[index].ishydrogen = ishydrogen;
}

AmoebaGKCavitationForce::NonbondedMethod AmoebaGKCavitationForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void AmoebaGKCavitationForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double AmoebaGKCavitationForce::getCutoffDistance() const {
    return cutoffDistance;
}

void AmoebaGKCavitationForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

void AmoebaGKCavitationForce::getParticleParameters(int index,  double& radius, double& gamma, double &vdw_alpha, double &charge,
				      bool& ishydrogen) const { 

    radius = particles[index].radius;
    gamma = particles[index].gamma;
    vdw_alpha = particles[index].vdw_alpha;
    charge = particles[index].charge;
    ishydrogen = particles[index].ishydrogen;
}

ForceImpl* AmoebaGKCavitationForce::createImpl() const {
    return new AmoebaGKCavitationForceImpl(*this);
}

void AmoebaGKCavitationForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaGKCavitationForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
