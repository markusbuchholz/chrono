// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Demo code illustrating the SCM semi-empirical model for deformable soil
// =============================================================================

#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChLinkMotorRotationAngle.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_irrlicht/ChIrrApp.h"

#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/terrain/SCMDeformableTerrain.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;
using namespace chrono::irrlicht;

using namespace irr;

bool output = false;
const std::string out_dir = GetChronoOutputPath() + "SCM_DEF_SOIL";

// Enable/disable adaptive mesh refinement
bool enable_adaptive_refinement = true;
double init_mesh_resolution = 0.1;
double min_mesh_resolution = 0.04;

// Enable/disable bulldozing effects
bool enable_bulldozing = true;

// Enable/disable moving patch feature
bool enable_moving_patch = false;

// If true, use provided callback to change soil properties based on location
bool var_params = true;

// Custom callback for setting location-dependent soil properties.
// Note that the (x,y) location is given in the terrain's reference plane. 
// Here, the vehicle moves in the terrain's negative y direction!
class MySoilParams : public vehicle::SCMDeformableTerrain::SoilParametersCallback {
  public:
    virtual void Set(double x, double y) override {
        if (y > 0) {
            m_Bekker_Kphi = 0.2e6;
            m_Bekker_Kc = 0;
            m_Bekker_n = 1.1;
            m_Mohr_cohesion = 0;
            m_Mohr_friction = 30;
            m_Janosi_shear = 0.01;
            m_elastic_K = 4e7;
            m_damping_R = 3e4;
        } else {
            m_Bekker_Kphi = 5301e3;
            m_Bekker_Kc = 102e3;
            m_Bekker_n = 0.793;
            m_Mohr_cohesion = 1.3e3;
            m_Mohr_friction = 31.1;
            m_Janosi_shear = 1.2e-2;
            m_elastic_K = 4e8;
            m_damping_R = 3e4;
        }
    }
};

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
    ChSystemSMC my_system;
    vehicle::SCMDeformableTerrain mterrain(&my_system);
    mterrain.Initialize(0, 6, 2, 60, 20, 2);
}
