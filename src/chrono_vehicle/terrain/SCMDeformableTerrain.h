#ifndef SCM_DEFORMABLE_GRID_TERRAIN_H
#define SCM_DEFORMABLE_GRID_TERRAIN_H

#include <set>
#include <string>
#include <unordered_map>
#include <ostream>

#include "chrono/assets/ChColorAsset.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/core/ChTimer.h"
#include "chrono/utils/ChConvexHull.h"


#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/ChSubsysDefs.h"
#include "chrono_vehicle/ChTerrain.h"

namespace chrono {
namespace vehicle {

// Utility Class and Functions for SCMDeformableGridTerrain
// =============================================================================
// ==========================faceVector.h=======================================
//==============================================================================
// faceVector object stores connection information for GridElement
// int a, b, c, d are the indexes for the vertices, and n_e is the normal vector
class faceVector {     
  public:
    int a;
    int b;
    int c;
    int d;
    int n_e;
};

// ============================================================================
// ===========================ChGridElement.h==================================
// ============================================================================
// Grid Element is one rectangular grid element
class ChGridElement{
  public:
    ChVector<> p1;  ///< first grid element vertex
    ChVector<> p2;  ///< second grid element vertex
    ChVector<> p3;  ///< third grid element vertex
    ChVector<> p4;  ///< fourth grid element vertex

    ChVector<> n;

  public:
    ChGridElement() : p1(VNULL), p2(VNULL), p3(VNULL), p4(VNULL), n(VNULL){}
    ChGridElement(const ChVector<>& mp1, const ChVector<>& mp2, const ChVector<>& mp3, const ChVector<>& mp4, const ChVector<>& mn): p1(mp1), p2(mp2), p3(mp3), p4(mp4), n(mn) {}
    ChGridElement(const ChGridElement& source);
    ~ChGridElement() {}

    /// Assignment operator: copy from another triangle
    ChGridElement& operator=(const ChGridElement& source);


    void GetBoundingBox(double& xmin,
                                double& xmax,
                                double& ymin,
                                double& ymax,
                                double& zmin,
                                double& zmax,
                                ChMatrix33<>* Rot = NULL);

    ChVector<> Baricenter();

    bool ContainVertex(ChVector<double> vertex);

};


// ======================================================================================
// ===========================ChSubGridMeshConnected.h==================================
// =======================================================================================
// sub mesh class
class ChSubGridMeshConnected{
  private:
    std::vector<ChGridElement> eleArr;
    std::vector<ChVector<double>> eleCenter;

    //std::vector<std::vector<int>> eleConnectedIdx;

    //mesh try -> store triangular mesh info
    std::vector<ChVector<>> mesh_vertices;
    std::vector<ChVector<int>> mesh_face;


//----------------------------------------- NEW
    std::vector<ChVector<>> vertices_vec;
    std::vector<std::vector<int>> neighbour_map_vec;
    std::vector<ChVector<int>> face_vec;
    std::vector<ChVector<>> colors_vec;


  public:
    void addGridElement(ChGridElement);
    void removeGridElement(int index);
    bool checkDataInteg();
    int getEleSize(){return eleArr.size();}
    int getCenterSize(){return eleCenter.size();}
    std::vector<ChGridElement> getEleArr(){return eleArr;}
    std::vector<ChVector<double>> getAllVertices();
    void getBoundingInfo();
    void Update(ChVector<double> new_vec, int idx);
    void Refine(ChVector<double> target_vertex);
    std::vector<ChVector<>> returnMeshVert();
    std::vector<ChVector<int>> returnMeshFace();
    void GetSubVisMesh(ChCoordsys<> plane);

    void UpdateColor(ChVector<float> new_color, int idx);

//--------------------------------------------------------
    void InitVerNeighMap();
    void InitAllVertices();
    void InitVerColor();

    std::vector<std::vector<int>> getAllNeigh_vec();
    std::vector<ChVector<double>> getAllVertices_vec();
    std::vector<ChVector<int>> getAllFaces();
    std::vector<ChVector<>> getAllColors_vec();
    

    double xmax;
    double xmin;
    double ymax;
    double ymin;
    double zmax;
    double zmin;

};

// =============================================================================
// ==========================ChGridMeshConnected.h=========================
// =============================================================================
// The entire mesh class
class ChGridMeshConnected{
  private:
    std::vector<ChSubGridMeshConnected> subArr;

  public:
    void initializeData(std::vector<ChGridElement> gird_ele, int sub_on_side);
    std::vector<ChSubGridMeshConnected> getSubGridData();
    std::vector<ChVector<double>> getAllVertices();
    void getBoundingInfo();
    void addSubGridData(ChSubGridMeshConnected subMesh);
    //std::shared_ptr<ChTriangleMeshShape> 
    void GetVisMesh(std::shared_ptr<ChTriangleMeshShape> trimesh,ChCoordsys<> plane, std::vector<int> active_sub_mesh,std::vector<int> vis_index);
    void Update(ChVector<double> new_vec,int idx ,int submesh_idx);
    void Refine(ChVector<double> target_vertex, int submesh_idx);
    void InitializeMeshVis(ChCoordsys<> plane);


    void UpdateColor(ChVector<float> new_color,int idx ,int submesh_idx);

    std::vector<ChSubGridMeshConnected>* getSubArrAddress(){
      return &subArr;
    }


    //------------------------
    void InitializeSubNeighMap();
    void InitSubAllVertices();
    void InitSubAllColors();
    


    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    

    std::vector<double> x_cut_Arr;
    std::vector<double> y_cut_Arr;



};


// =============================================================================
// ==========================GridMeshLoader.h===================================
// =============================================================================
// Helper class - Load a rectangular mesh into the program
// The return is a vector of GridElement
class GridMeshLoader {     
  private:
    std::vector<ChVector<>> vertices;
    std::vector<faceVector> faces;
    std::vector<ChVector<>> normals;
    std::vector<ChGridElement> gridEle;
    bool loadFromGridWaveFront(std::string path);

  public:             
    
    void addVertice(ChVector<>& a);
    std::vector<ChGridElement> loadObj(std::string path);
};


/// @addtogroup vehicle_terrain
class SCMDeformableSoilGrid;
class CH_VEHICLE_API SCMDeformableTerrain : public ChTerrain {
  private:
      std::shared_ptr<ChGridMeshConnected> Grid;
      std::shared_ptr<SCMDeformableSoilGrid> m_ground;

  public:
    enum DataPlotType {
        PLOT_NONE,
        PLOT_LEVEL,
        PLOT_LEVEL_INITIAL,
        PLOT_SINKAGE,
        PLOT_SINKAGE_ELASTIC,
        PLOT_SINKAGE_PLASTIC,
        PLOT_STEP_PLASTIC_FLOW,
        PLOT_PRESSURE,
        PLOT_PRESSURE_YELD,
        PLOT_SHEAR,
        PLOT_K_JANOSI,
        PLOT_IS_TOUCHED,
        PLOT_ISLAND_ID,
        PLOT_MASSREMAINDER
    };

    TerrainForce GetContactForce(std::shared_ptr<ChBody> body) const;

    void SetPlane(ChCoordsys<> mplane);

    void SetPlotType(DataPlotType mplot, double mmin, double mmax);

      SCMDeformableTerrain(ChSystem* system,bool visualization_mesh = true);

      ~SCMDeformableTerrain(){}

      void Initialize(double height, double sizeX, double sizeZ, int divX, int divZ, int sub_per_side=2, ChCoordsys<> plane=ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(-CH_C_PI_2)));

      void Initialize(const std::string& mesh_file, int sub_per_side=2, ChCoordsys<> plane=ChCoordsys<>(ChVector<>(0, 0, 0), Q_from_AngX(-CH_C_PI_2)));

      void SetSoilParameters(
          double Bekker_Kphi,    // Kphi, frictional modulus in Bekker model
          double Bekker_Kc,      // Kc, cohesive modulus in Bekker model
          double Bekker_n,       // n, exponent of sinkage in Bekker model (usually 0.6...1.8)
          double Mohr_cohesion,  // Cohesion in, Pa, for shear failure
          double Mohr_friction,  // Friction angle (in degrees!), for shear failure
          double Janosi_shear,   // J , shear parameter, in meters, in Janosi-Hanamoto formula (usually few mm or cm)
          double elastic_K,      // elastic stiffness K (must be > Kphi; very high values gives the original SCM model)
          double damping_R  // vertical damping R, per unit area (vertical speed proportional, it is zero in original SCM model)
      );

      void SetBulldozingParameters(
          double mbulldozing_erosion_angle,       ///< angle of erosion of the displaced material (in degrees!)
          double mbulldozing_flow_factor,         ///< growth of lateral volume respect to pressed volume
          int mbulldozing_erosion_n_iterations,   ///< number of erosion refinements per timestep
          int mbulldozing_erosion_n_propagations  ///< number of concentric vertex selections subject to erosion
      );

      void SetBulldozingFlow(bool mb);
      bool GetBulldozingFlow() const;


      void SetAutomaticRefinement(bool mr);
      bool GetAutomaticRefinement();


      // this function might not be needed anymore
      void SetAutomaticRefinementResolution(double lv);
      // this function might not be needed anymore
      int GetAutomaticRefinementResolution();

      void SetTestHighOffset(double mr);
      double GetTestHighOffset();
      void AddMovingPatch(std::shared_ptr<ChBody> body,
                                             const ChVector<>& point_on_body,
                                             double dimX,
                                             double dimY);

      int returnVerticesSize();
      std::vector<ChVector<double>> returnVertices();

      void PrintStepStatistics(std::ostream& os) const;

      /// Get the underlying triangular mesh.
      const std::shared_ptr<ChTriangleMeshShape> GetMesh() const;

      void SetColor(ChColor color);




    /// Class to be used as a callback interface for location-dependent soil parameters.
    /// A derived class must implement Set() and set **all** soil parameters (no defaults are provided).
    class CH_VEHICLE_API SoilParametersCallback {
      public:
        virtual ~SoilParametersCallback() {}

        /// Set the soil properties at a given (x,y) location.
        /// Attention: the location is assumed to be provided in the (x,y) plane of the patch reference plane!
        /// An implementation of this method in a derived class must set all soil parameters.
        virtual void Set(double x, double y) = 0;

        double m_Bekker_Kphi;    ///< Kphi, frictional modulus in Bekker model
        double m_Bekker_Kc;      ///< Kc, cohesive modulus in Bekker model
        double m_Bekker_n;       ///< n, exponent of sinkage in Bekker model (usually 0.6...1.8)
        double m_Mohr_cohesion;  ///< Cohesion in, Pa, for shear failure
        double m_Mohr_friction;  ///< Friction angle (in degrees!), for shear failure
        double m_Janosi_shear;   ///< J , shear parameter, in meters, in Janosi-Hanamoto formula (usually few mm or cm)
        double m_elastic_K;      ///< elastic stiffness K, per unit area, [Pa/m] (must be larger than Kphi)
        double m_damping_R;      ///< vertical damping R, per unit area [Pa s/m] (proportional to vertical speed)
    };
      void RegisterSoilParametersCallback(std::shared_ptr<SoilParametersCallback> cb);


};


class CH_VEHICLE_API SCMDeformableSoilGrid : public ChLoadContainer {
  public:
       struct MovingPatchInfo {
        std::shared_ptr<ChBody> m_body;  // tracked body
        ChVector<> m_point;              // patch center, relative to body
        ChVector2<> m_dim;               // patch dimensions (X,Y)
        ChVector2<> m_min;               // current patch AABB (min x,y)
        ChVector2<> m_max;               // current patch AABB (max x,y)
    };


    SCMDeformableSoilGrid(ChSystem* system, bool visualization_mesh);
    ~SCMDeformableSoilGrid() {}

    /// Initialize the terrain system with a grid mesh file
    void Initialize(std::shared_ptr<ChGridMeshConnected> Grid);
    void GetMesh();


  private:

    // Get the terrain height below the specified location.
    double GetHeight(const ChVector<>& loc) const;

    void UpdateFixedPatch();

    // Updates the forces and the geometry, at the beginning of each timestep
    virtual void Setup() override {
        this->ComputeInternalForces();
        ChLoadContainer::Update(ChTime, true);
    }

    // Updates the forces and the geometry
    virtual void Update(double mytime, bool update_assets = true) override {
        ChTime = mytime;
    }

    // Reset the list of forces, and fills it with forces from a soil contact model.
    // This is called automatically during timestepping (only at the beginning of
    // each IntLoadResidual_F() for performance reason, not at each Update() that might be overkill).
    void ComputeInternalForces();

    // Override the ChLoadContainer method for computing the generalized force F term:
    virtual void IntLoadResidual_F(const unsigned int off,  ///< offset in R residual
                                   ChVectorDynamic<>& R,    ///< result: the R residual, R += c*F
                                   const double c           ///< a scaling factor
                                   ) override {
        // Overloading base class, that takes all F vectors from the list of forces and put all them in R
        ChLoadContainer::IntLoadResidual_F(off, R, c);
    }

    // This is called after Initialize(), it pre-computes aux.topology
    // data structures for the mesh, aux. material data, etc.
    void SetupAuxData();


    std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> y_cut, 
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches);

    std::shared_ptr<ChColorAsset> m_color;
    std::shared_ptr<ChGridMeshConnected> m_grid_shape;
    double m_height;

    std::vector<int> active_sub_mesh;

    std::vector<ChVector<>> vertices;

    std::vector<ChVector<>> p_vertices_initial;
    std::vector<ChVector<>> p_speeds;
    std::vector<double> p_level;
    std::vector<double> p_level_initial;
    std::vector<double> p_hit_level;
    std::vector<double> p_sinkage;
    std::vector<double> p_sinkage_plastic;
    std::vector<double> p_sinkage_elastic;
    std::vector<double> p_step_plastic_flow;
    std::vector<double> p_kshear;  // Janosi-Hanamoto shear accumulator
    std::vector<double> p_area;
    std::vector<double> p_sigma;
    std::vector<double> p_sigma_yeld;
    std::vector<double> p_tau;
    std::vector<double> p_massremainder;
    std::vector<int> p_id_island;
    std::vector<bool> p_erosion;

    double m_Bekker_Kphi;
    double m_Bekker_Kc;
    double m_Bekker_n;
    double m_Mohr_cohesion;
    double m_Mohr_friction;
    double m_Janosi_shear;
    double m_elastic_K;
    double m_damping_R;

    int plot_type;
    double plot_v_min;
    double plot_v_max;

    double area_x;
    double area_y;

    ChCoordsys<> plane;
    //PatchType m_type;

    // aux. topology data
    //std::vector<std::set<int>> connected_vertexes;
    std::shared_ptr<ChTriangleMeshShape> m_trimesh_shape;
    //std::vector<std::array<int, 4>> tri_map;

    bool do_bulldozing;
    double bulldozing_flow_factor;
    double bulldozing_erosion_angle;
    int bulldozing_erosion_n_iterations;
    int bulldozing_erosion_n_propagations;

    bool do_refinement;
    double refinement_level;

    double test_high_offset;
    double test_low_offset;

    double last_t;  // for optimization




    //this is the parmeter for findActiveIdx
    double dx_findactive; //submesh x dim
    double dy_findactive; //submesh y dim
    double dside_findactive; //division per side
    double tot_findactive; //total number of submeshes



    // this vector stores the index cut data for mesh visualization
    std::vector<int> vis_index;




    // Moving patch parameters
 
    std::vector<MovingPatchInfo> m_patches;  // set of active moving patches
    bool m_moving_patch;                     // moving patch feature enabled?

    // Callback object for position-dependent soil properties
    std::shared_ptr<SCMDeformableTerrain::SoilParametersCallback> m_soil_fun;

    // Timers and counters
    ChTimer<double> m_timer_calc_areas;
    ChTimer<double> m_timer_ray_casting;
    ChTimer<double> m_timer_update;
    ChTimer<double> m_timer_refinement;
    ChTimer<double> m_timer_bulldozing;
    ChTimer<double> m_timer_visualization;
    ChTimer<double> m_timer_total;
    ChTimer<double> m_timer_vertsetup;
    ChTimer<double> m_timer_loadlist;
    
    
    size_t m_num_vertices;
    size_t m_num_faces;
    size_t m_num_ray_casts;
    size_t m_num_marked_faces;

    std::unordered_map<ChContactable*, TerrainForce> m_contact_forces;

    std::shared_ptr<ChGridMeshConnected> GetGrid(){return m_grid_shape;}

    friend class SCMDeformableTerrain;
};

/// @} vehicle_terrain


}  // end namespace vehicle
}  // end namespace chrono

#endif
