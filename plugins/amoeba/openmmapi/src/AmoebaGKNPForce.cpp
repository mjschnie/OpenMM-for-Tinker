/* -------------------------------------------------------------------------- *
 *                             OpenMM-GKNP                                  *
 * -------------------------------------------------------------------------- */

#include <iostream>
#include "openmm/AmoebaGKNPForce.h"
#include "openmm/internal/AmoebaGKNPForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace OpenMM;
using namespace std;

AmoebaGKNPForce::AmoebaGKNPForce() : nonbondedMethod(NoCutoff), cutoffDistance(1.0), solvent_radius(SOLVENT_RADIUS) {
}

int AmoebaGKNPForce::addParticle(double radius, double gamma, double vdw_alpha, double charge, bool ishydrogen){
  ParticleInfo particle(radius, gamma, vdw_alpha, charge, ishydrogen);
  particles.push_back(particle);
  return particles.size()-1;
}

void AmoebaGKNPForce::setParticleParameters(int index, double radius, double gamma, double vdw_alpha, double charge, bool ishydrogen){
  ASSERT_VALID_INDEX(index, particles);
  particles[index].radius = radius;
  particles[index].radius;
  particles[index].gamma = gamma;
  particles[index].vdw_alpha = vdw_alpha;
  particles[index].charge = charge;
  particles[index].ishydrogen = ishydrogen;
}

AmoebaGKNPForce::NonbondedMethod AmoebaGKNPForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void AmoebaGKNPForce::setNonbondedMethod(NonbondedMethod method) {
    nonbondedMethod = method;
}

double AmoebaGKNPForce::getCutoffDistance() const {
    return cutoffDistance;
}

void AmoebaGKNPForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

void AmoebaGKNPForce::getParticleParameters(int index,  double& radius, double& gamma, double &vdw_alpha, double &charge,
				      bool& ishydrogen) const { 

    ASSERT_VALID_INDEX(index, particles);
    radius = particles[index].radius;
    gamma = particles[index].gamma;
    vdw_alpha = particles[index].vdw_alpha;
    charge = particles[index].charge;
    ishydrogen = particles[index].ishydrogen;
}

ForceImpl* AmoebaGKNPForce::createImpl() const {
    return new AmoebaGKNPForceImpl(*this);
}

void AmoebaGKNPForce::updateParametersInContext(Context& context) {
    dynamic_cast<AmoebaGKNPForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
