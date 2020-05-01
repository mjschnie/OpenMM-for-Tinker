#ifndef OPENMM_AMOEBAREFERENCEGKCAVITATIONFORCE_H
#define OPENMM_AMOEBAREFERENCEGKCAVITATIONFORCE_H

#include "openmm/Vec3.h"
#include <string>
#include <vector>
#include <openmm/Platform.h>
#include <openmm/System.h>
#include <openmm/AmoebaGKCavitationForce.h>
#include "gaussvol.h"


namespace OpenMM {
    class AmoebaReferenceGKCavitationForce{
    public:
        AmoebaReferenceGKCavitationForce(){
            gvol = 0;
        };
        ~AmoebaReferenceGKCavitationForce(){
            if(gvol) delete gvol;
        };
        double calculateForceAndEnergy(vector<RealVec>& pos, vector<RealVec>& force, int numParticles,
                vector<int> ishydrogen, vector<RealOpenMM> radii_large, vector<RealOpenMM> radii_vdw, vector<RealOpenMM> gammas,
                double roffset, vector<RealVec> vol_force, vector<RealOpenMM> vol_dv, vector<RealOpenMM> free_volume,
                vector<RealOpenMM> self_volume);
    private:
        GaussVol *gvol; // gaussvol instance
    };
}

#endif //OPENMM_AMOEBAREFERENCEGKCAVITATIONFORCE_H
