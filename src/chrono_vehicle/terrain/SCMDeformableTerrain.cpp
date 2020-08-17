
#include "chrono_vehicle/terrain/SCMDeformableTerrain.h"
#include <math.h> 
#include <queue>


namespace chrono {
namespace vehicle {

// Utility Classes for SCDeformableGridTerrain
// ====================================================================
// ==========================ChGridElement.cpp=========================
// ====================================================================
ChGridElement::ChGridElement(const ChGridElement& source) {
    p1 = source.p1; // the first vertex of a rectangular element
    p2 = source.p2; // the second vertex of a rectangular element
    p3 = source.p3; // the third vertex of a rectangular element
    p4 = source.p4; // the fourth vertex of a rectangular element
    n = source.n;   // the normal of a rectangular element
}

ChGridElement& ChGridElement::operator=(const ChGridElement& source) {
    if (&source == this)
        return *this;
    p1 = source.p1;
    p2 = source.p2;
    p3 = source.p3;
    p4 = source.p4;
    n = source.n;
    return *this;
}

// Return the bounding box of one specific rectangular element
void ChGridElement::GetBoundingBox(double& xmin,
                                double& xmax,
                                double& ymin,
                                double& ymax,
                                double& zmin,
                                double& zmax,
                                ChMatrix33<>* Rot){
    if (Rot == NULL) {
        xmin = ChMin(ChMin(ChMin(p1.x(), p2.x()), p3.x()), p4.x());
        ymin = ChMin(ChMin(ChMin(p1.y(), p2.y()), p3.y()), p4.y());
        zmin = ChMin(ChMin(ChMin(p1.z(), p2.z()), p3.z()), p4.z());
        xmax = ChMax(ChMax(ChMax(p1.x(), p2.x()), p3.x()), p4.x());
        ymax = ChMax(ChMax(ChMax(p1.y(), p2.y()), p3.y()), p4.y());
        zmax = ChMax(ChMax(ChMax(p1.z(), p2.z()), p3.z()), p4.z());
    } else {
        ChVector<> trp1 = Rot->transpose() * p1;
        ChVector<> trp2 = Rot->transpose() * p2;
        ChVector<> trp3 = Rot->transpose() * p3;
        ChVector<> trp4 = Rot->transpose() * p4;
        xmin = ChMin(ChMin(ChMin(p1.x(), p2.x()), p3.x()), p4.x());
        ymin = ChMin(ChMin(ChMin(p1.y(), p2.y()), p3.y()), p4.y());
        zmin = ChMin(ChMin(ChMin(p1.z(), p2.z()), p3.z()), p4.z());
        xmax = ChMax(ChMax(ChMax(p1.x(), p2.x()), p3.x()), p4.x());
        ymax = ChMax(ChMax(ChMax(p1.y(), p2.y()), p3.y()), p4.y());
        zmax = ChMax(ChMax(ChMax(p1.z(), p2.z()), p3.z()), p4.z());
    }
}

// return the center of the grid element
ChVector<> ChGridElement::Baricenter(){
    ChVector<> mb;
    mb.x() = (p1.x() + p2.x() + p3.x() + p4.x()) / 4.;
    mb.y() = (p1.y() + p2.y() + p3.y() + p4.y()) / 4.;
    mb.z() = (p1.z() + p2.z() + p3.z() + p4.z()) / 4.;
    return mb;
}

// =============================================================================
// ==========================ChSubGridMeshConnected.cpp=========================
// =============================================================================

// helper function check repeat, check whether a ChVector<doouble> is already in the array
bool checkRepeat(ChVector<double> a, std::vector<ChVector<double>> arr);

// add a new rectangular element in the sub mesh
void ChSubGridMeshConnected::addGridElement(ChGridElement temp){
    eleArr.push_back(temp);
    ChVector<double> buffer = temp.Baricenter();
    eleCenter.push_back(buffer);
}

// remove a rectangular element in the sub mesh by index
void ChSubGridMeshConnected::removeGridElement(int index){
    eleArr.erase(eleArr.begin()+index);
    eleCenter.erase(eleCenter.begin()+index);
}

// determine whether the center data is complete by
// checking whether the number of baricenter = the number of rectangular elements in sub mesh
bool ChSubGridMeshConnected::checkDataInteg(){
    int eleArrSize = eleArr.size();
    int eleCenterSize = eleCenter.size();

    if(eleArrSize==eleCenterSize){
        return true;
    }else{
        return false;
    }
}

// return all vertices in a submesh
std::vector<ChVector<double>> ChSubGridMeshConnected::getAllVertices(){
    std::vector<ChVector<double>> returnArr;
    // index indicator, so we can get rectangular connection infomation after
    // getAllVertices() is called
    // int index_ind = 0;  

    for(int i = 0; i<eleArr.size();i++){

        if(checkRepeat(eleArr[i].p1, returnArr)==false){
            returnArr.push_back(eleArr[i].p1);
        }

        if(checkRepeat(eleArr[i].p2, returnArr)==false){
            returnArr.push_back(eleArr[i].p2);
        }

        if(checkRepeat(eleArr[i].p3, returnArr)==false){
            returnArr.push_back(eleArr[i].p3);
        }

        if(checkRepeat(eleArr[i].p4, returnArr)==false){
            returnArr.push_back(eleArr[i].p4);
        }
    }
    return returnArr;
}

// return the bounding box info of a specific submesh
void ChSubGridMeshConnected::getBoundingInfo(){
    double xminBuff = 999.;
    double xmaxBuff = -999.;
    double yminBuff = 999.;
    double ymaxBuff = -999.;
    double zminBuff = 999.;
    double zmaxBuff = -999.;
    for(int i = 0; i<eleArr.size();i++){
        double xmintemp;
        double xmaxtemp;
        double ymintemp;
        double ymaxtemp;
        double zmintemp;
        double zmaxtemp;

        eleArr[i].GetBoundingBox(xmintemp,xmaxtemp,ymintemp,ymaxtemp,zmintemp,zmaxtemp);

        if(xmintemp < xminBuff){xminBuff = xmintemp;}
        if(xmaxtemp > xmaxBuff){xmaxBuff = xmaxtemp;}
        if(ymintemp < yminBuff){yminBuff = ymintemp;}
        if(ymaxtemp > ymaxBuff){ymaxBuff = ymaxtemp;}
        if(zmintemp < zminBuff){zminBuff = zmintemp;}
        if(zmaxtemp > zmaxBuff){zmaxBuff = zmaxtemp;}
    }

    // update max and min info
    xmax = xmaxBuff;
    xmin = xminBuff;
    ymax = ymaxBuff;
    ymin = yminBuff;
    zmax = ymaxBuff;
    zmin = zminBuff;
}

bool checkRepeat(ChVector<double> a, std::vector<ChVector<double>> arr){
    for(int i = 0; i<arr.size();i++){
        if(arr[i] == a){
            return true;
        }
    }

    return false;
}

// =============================================================================
// ==========================ChGridMeshConnected.cpp============================
//==============================================================================
void ChGridMeshConnected::initializeData(std::vector<ChGridElement> grid_ele, int sub_on_side){

    for (int i = 0; i<grid_ele.size();i++){
        std::cout<<"ele: "<<i<<std::endl;
        std::cout<<"p1: "<<grid_ele[i].p1<<std::endl;
        std::cout<<"p2: "<<grid_ele[i].p2<<std::endl;
        std::cout<<"p3: "<<grid_ele[i].p3<<std::endl;
        std::cout<<"p4: "<<grid_ele[i].p4<<std::endl;
        std::cout<<"=================="<<std::endl;
    }




    double x_max=-999.;
    double x_min=999.;
    double y_max=-999.;
    double y_min=999.;
    double z_max=-999.;
    double z_min=999.;
    for(int i = 0; i<grid_ele.size();i++){
        double x_max_buffer;
        double x_min_buffer;
        double y_max_buffer;
        double y_min_buffer;
        double z_max_buffer;
        double z_min_buffer;
        grid_ele[i].GetBoundingBox(x_min_buffer,x_max_buffer,y_min_buffer,y_max_buffer,z_min_buffer,z_max_buffer);
    
    if(x_min_buffer<x_min){
        x_min = x_min_buffer;
    }

    if(x_max_buffer>x_max){
        x_max = x_max_buffer;
    }

    if(y_min_buffer<y_min){
        y_min = y_min_buffer;
    }

    if(y_max_buffer>y_max){
        y_max = y_max_buffer;
    }

    if(z_min_buffer<z_min){
        z_min = z_min_buffer;
    }

    if(z_max_buffer>z_max){
        z_max = z_max_buffer;
    }
    }

    double x_tot_dis = x_max - x_min;
    double z_tot_dis = z_max - z_min;

    std::cout<<"x_tot_dis: "<<x_tot_dis<<std::endl;
    std::cout<<"z_tot_dis: "<<z_tot_dis<<std::endl;

    std::vector<double> x_cut;
    std::vector<double> z_cut;

    for(int i = 0; i<sub_on_side; i++){
        x_cut.push_back(x_min+(x_tot_dis/sub_on_side)*(i+1));
        z_cut.push_back(z_min+(z_tot_dis/sub_on_side)*(i+1));
    }

    // add GridElements to SubMesh and construct a connected Grid Mesh structure
    for(int i = 0; i<x_cut.size();i++){
        for(int j = 0; j<z_cut.size();j++){
            ChSubGridMeshConnected subTemp;

            for(int k = 0;k<grid_ele.size();k++){
                if(grid_ele[k].Baricenter().x()<x_cut[i] && grid_ele[k].Baricenter().z()<z_cut[j]){
                    subTemp.addGridElement(grid_ele[k]);
                    grid_ele.erase(grid_ele.begin()+k);
                    k=k-1;
                }
               
            }
            addSubGridData(subTemp);
        }
    }

    
}

void ChGridMeshConnected::addSubGridData(ChSubGridMeshConnected subMesh){
    subArr.push_back(subMesh);
}


std::vector<ChSubGridMeshConnected> ChGridMeshConnected::getSubGridData(){
    // return a vector of sub mesh
    return subArr;
}

std::vector<ChVector<double>> ChGridMeshConnected::getAllVertices(){
    std::vector<ChVector<double>> returnArr;

    // might need to check repeat on edges?
    for(int i = 0; i<subArr.size();i++){
        std::vector<ChVector<double>> returnSubArr;
        returnSubArr = subArr[i].getAllVertices();
        returnArr.insert(returnArr.end(), returnSubArr.begin(), returnSubArr.end());
    }
    return returnArr;
}

void ChGridMeshConnected::getBoundingInfo(){
    double xminBuff = 999.;
    double xmaxBuff = -999;
    double yminBuff = 999.;
    double ymaxBuff = -999;
    double zminBuff = 999;
    double zmaxBuff = -999;
    for(int i = 0; i<subArr.size();i++){
        double xmintemp = subArr[i].xmin;
        double xmaxtemp = subArr[i].xmax;
        double ymintemp = subArr[i].ymin;
        double ymaxtemp = subArr[i].ymax;
        double zmintemp = subArr[i].zmin;
        double zmaxtemp = subArr[i].zmax;

        if(xmintemp < xminBuff){xminBuff = xmintemp;}
        if(xmaxtemp > xmaxBuff){xmaxBuff = xmaxtemp;}
        if(ymintemp < yminBuff){yminBuff = ymintemp;}
        if(ymaxtemp > ymaxBuff){ymaxBuff = ymaxtemp;}
        if(zmintemp < zminBuff){zminBuff = zmintemp;}
        if(zmaxtemp > zmaxBuff){zmaxBuff = zmaxtemp;}
    }

    xmax = xmaxBuff;
    xmin = xminBuff;
    ymax = ymaxBuff;
    ymin = yminBuff;
    zmax = ymaxBuff;
    zmin = zminBuff;
}


// =============================================================================
// ==========================GridMeshLoader.cpp===================================
// ==============================================================================
// Helper class to load a grid mesh
std::vector<std::string> splitHelper(const std::string& s, char seperator);
std::vector<std::string> split(std::string stringToBeSplitted, std::string delimeter);

void GridMeshLoader::addVertice(ChVector<>& a){
    vertices.push_back(a);
}

bool GridMeshLoader::loadFromGridWaveFront(std::string path){
    std::ifstream file(path);

    if (!file.is_open())
	    return false;


    std::string curline;
    std::string::size_type sz;
    while(std::getline(file, curline)){
        
        if (curline.rfind("v", 0)==0) {
            std::vector<std::string> buffer = splitHelper(curline,' ');
            
            double a1 = std::stod(buffer[1],&sz);
            double a2 = std::stod(buffer[2],&sz);
            double a3 = std::stod(buffer[3],&sz);

            vertices.push_back(ChVector<>({a1, a2, a3}));
        }

        if(curline.rfind("f",0)==0){
            std::vector<std::string> buffer = splitHelper(curline,' ');

            int f1;
            int f2;
            int f3;
            int f4;
            int vn_index;

            for(int i = 1;i<buffer.size();i++){
                
                std::vector<std::string> buffer_buffer = split(buffer[i],"//");
                if(i==1){
                    f1 = std::stod(buffer_buffer[0],&sz);
                    
                    vn_index = std::stod(buffer_buffer[1],&sz);
                }

                if(i==2){
                    f2 = std::stod(buffer_buffer[0],&sz);
                }

                if(i==3){
                    f3 = std::stod(buffer_buffer[0],&sz);
                }

                if(i==4){
                    f4 = std::stod(buffer_buffer[0],&sz);
                }
            }

            faceVector temp;
            temp.a = f1;
            temp.b = f2;
            temp.c = f3;
            temp.d = f4;
            temp.n_e = vn_index;

            faces.push_back(temp);
        }

        if(curline.rfind("vn",0)==0){
            std::vector<std::string> buffer = splitHelper(curline,' ');
            double a1 = std::stod(buffer[1],&sz);
            double a2 = std::stod(buffer[2],&sz);
            double a3 = std::stod(buffer[3],&sz);
            normals.push_back(ChVector<>({a1, a2, a3}));
        }
    }
    
    for(int i = 0; i<faces.size();i++){
        gridEle.push_back(ChGridElement(vertices[faces[i].a-1], vertices[faces[i].b-1],
            vertices[faces[i].c-1], vertices[faces[i].d-1], normals[faces[i].n_e-1]));
    }

    file.close();

    return true;
}

std::vector<ChGridElement> GridMeshLoader::loadObj(std::string path){
    loadFromGridWaveFront(path);
    return gridEle;
    
}

// Helper class to split a string into 2 by char seperator, this is for mesh loading
std::vector<std::string> splitHelper(const std::string& s, char seperator)
{
   std::vector<std::string> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(substring);

        prev_pos = ++pos;
    }

    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word

    return output;
}

// Helper class to split a string into 2 by string seperator, this is for mesh loading
std::vector<std::string> split(std::string stringToBeSplitted, std::string delimeter)
{
     std::vector<std::string> splittedString;
     int startIndex = 0;
     int  endIndex = 0;
     while( (endIndex = stringToBeSplitted.find(delimeter, startIndex)) < stringToBeSplitted.size() )
    {
       std::string val = stringToBeSplitted.substr(startIndex, endIndex - startIndex);
       splittedString.push_back(val);
       startIndex = endIndex + delimeter.size();
     }
     if(startIndex < stringToBeSplitted.size())
     {
       std::string val = stringToBeSplitted.substr(startIndex);
       splittedString.push_back(val);
     }
     return splittedString;
}


// Helper class to find all active sub meshes
// This function returns a vector of int values 
// which are the indexes of active submesh
std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> z_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<std::vector<chrono::vehicle::SCMDeformableSoilGrid::MovingPatchInfo>
>patches);

// Check whether a int value already exists in cut_cross vector
bool CheckIdxRepeat(int target,std::vector<int> cut_cross);

SCMDeformableTerrain::SCMDeformableTerrain(ChSystem* system, bool visualization_mesh ) {
    m_ground = chrono_types::make_shared<SCMDeformableSoilGrid>(system, visualization_mesh);
    system->Add(m_ground);
    Grid = std::make_shared<ChGridMeshConnected>();
}

// Initialize a grid mesh by input values
void SCMDeformableTerrain::Initialize(double height, double sizeX, double sizeZ, int divX, int divZ, int sub_per_side)
{
    std::vector<double> x_cut;
    std::vector<double> z_cut;

    double dx = sizeX / divX;
    double dz = sizeZ / divZ;

    std::vector<ChGridElement> temp_store_buffer;

    for(double i = 0; i<divX; i++){
        for(double j = 0; j<divZ; j++){
            ChVector<double> p1;
            p1.Set(0+i*dx, height, 0+j*dz);
            ChVector<double> p2;
            p2.Set(0+(i+1)*dx, height, 0+j*dz);
            ChVector<double> p3;
            p3.Set(0+i*dx, height, 0+(j+1)*dz);
            ChVector<double> p4;
            p4.Set(0+(i+1)*dx, height, 0+(j+1)*dz);
            ChVector<double> n;
            n.Set(0, 1, 0);
            temp_store_buffer.push_back(ChGridElement(p1, p2, p3, p4,n));
        }
    }
    //Grid = new ChGridMeshConnected();
    Grid->initializeData(temp_store_buffer, sub_per_side);
    std::cout<<"test point"<<std::endl;
    m_ground->Initialize(Grid);
}




void SCMDeformableTerrain::Initialize(const std::string& mesh_file, int sub_per_side)
{
    GridMeshLoader test;
    std::vector<ChGridElement> result = test.loadObj(mesh_file);
    Grid->initializeData(result, sub_per_side);
    m_ground->Initialize(Grid);
}




void SCMDeformableTerrain::SetSoilParameters(
    double Bekker_Kphi,    // Kphi, frictional modulus in Bekker model
    double Bekker_Kc,      // Kc, cohesive modulus in Bekker model
    double Bekker_n,       // n, exponent of sinkage in Bekker model (usually 0.6...1.8)
    double Mohr_cohesion,  // Cohesion in, Pa, for shear failure
    double Mohr_friction,  // Friction angle (in degrees!), for shear failure
    double Janosi_shear,   // J , shear parameter, in meters, in Janosi-Hanamoto formula (usually few mm or cm)
    double elastic_K,      // elastic stiffness K (must be > Kphi; very high values gives the original SCM model)
    double damping_R  // vertical damping R, per unit area (vertical speed proportional, it is zero in original SCM model)
) {
    m_ground->m_Bekker_Kphi = Bekker_Kphi;
    m_ground->m_Bekker_Kc = Bekker_Kc;
    m_ground->m_Bekker_n = Bekker_n;
    m_ground->m_Mohr_cohesion = Mohr_cohesion;
    m_ground->m_Mohr_friction = Mohr_friction;
    m_ground->m_Janosi_shear = Janosi_shear;
    m_ground->m_elastic_K = ChMax(elastic_K, Bekker_Kphi);
    m_ground->m_damping_R = damping_R;
}




void SCMDeformableTerrain::SetBulldozingParameters(
    double mbulldozing_erosion_angle,       ///< angle of erosion of the displaced material (in degrees!)
    double mbulldozing_flow_factor,         ///< growth of lateral volume respect to pressed volume
    int mbulldozing_erosion_n_iterations,   ///< number of erosion refinements per timestep
    int mbulldozing_erosion_n_propagations  ///< number of concentric vertex selections subject to erosion
) {
    m_ground->bulldozing_erosion_angle = mbulldozing_erosion_angle;
    m_ground->bulldozing_flow_factor = mbulldozing_flow_factor;
    m_ground->bulldozing_erosion_n_iterations = mbulldozing_erosion_n_iterations;
    m_ground->bulldozing_erosion_n_propagations = mbulldozing_erosion_n_propagations;
}

void SCMDeformableTerrain::SetAutomaticRefinement(bool mr) {
    m_ground->do_refinement = mr;
}

bool SCMDeformableTerrain::GetAutomaticRefinement(){
    return m_ground->do_refinement;
}

void SCMDeformableTerrain::SetAutomaticRefinementResolution(double lv) {
    m_ground->refinement_level = lv;
}

int SCMDeformableTerrain::GetAutomaticRefinementResolution(){
    return m_ground->refinement_level;
}

void SCMDeformableTerrain::SetTestHighOffset(double mr) {
    m_ground->test_high_offset = mr;
}

double SCMDeformableTerrain::GetTestHighOffset(){
    return m_ground->test_high_offset;
}

void SCMDeformableTerrain::SetBulldozingFlow(bool mb) {
    m_ground->do_bulldozing = mb;
}

bool SCMDeformableTerrain::GetBulldozingFlow() const {
    return m_ground->do_bulldozing;
}

const std::shared_ptr<ChTriangleMeshShape> SCMDeformableTerrain::GetMesh() const {
    return m_ground->m_trimesh_shape;
}

void SCMDeformableTerrain::AddMovingPatch(std::shared_ptr<ChBody> body,
                                             const ChVector<>& point_on_body,
                                             double dimX,
                                             double dimY) {
    SCMDeformableSoilGrid::MovingPatchInfo pinfo;
    pinfo.m_body = body;
    pinfo.m_point = point_on_body;
    pinfo.m_dim = ChVector2<>(dimX, dimY);

    m_ground->m_patches.push_back(pinfo);

    // Moving patch monitoring is now enabled
    m_ground->m_moving_patch = true;
}

// MIGHT NEED TO BE REMOVED?
void SCMDeformableTerrain::SetPlotType(DataPlotType mplot, double mmin, double mmax) {
    m_ground->plot_type = mplot;
    m_ground->plot_v_min = mmin;
    m_ground->plot_v_max = mmax;
}

int SCMDeformableTerrain::returnVerticesSize(){
    return Grid->getAllVertices().size();
}

std::vector<ChVector<double>> SCMDeformableTerrain::returnVertices(){
    return Grid->getAllVertices();
}

// Set the plane reference.
void SCMDeformableTerrain::SetPlane(ChCoordsys<> mplane) {
    m_ground->plane = mplane;
}

void SCMDeformableTerrain::PrintStepStatistics(std::ostream& os) const {
    os << " Timers:" << std::endl;
    os << "   Calculate areas:         " << m_ground->m_timer_calc_areas() << std::endl;
    os << "   Ray casting:             " << m_ground->m_timer_ray_casting() << std::endl;
    if (m_ground->do_refinement)
        os << "   Refinements:             " << m_ground->m_timer_refinement() << std::endl;
    if (m_ground->do_bulldozing)
        os << "   Bulldozing:              " << m_ground->m_timer_bulldozing() << std::endl;
    os << "   Visualization:           " << m_ground->m_timer_visualization() << std::endl;

    os << " Counters:" << std::endl;
    os << "   Number vertices:         " << m_ground->m_num_vertices << std::endl;
    os << "   Number ray-casts:        " << m_ground->m_num_ray_casts << std::endl;
    os << "   Number faces:            " << m_ground->m_num_faces << std::endl;
    if (m_ground->do_refinement)
        os << "   Number faces refinement: " << m_ground->m_num_marked_faces << std::endl;
}



TerrainForce SCMDeformableTerrain::GetContactForce(std::shared_ptr<ChBody> body) const {
    auto itr = m_ground->m_contact_forces.find(body.get());
    if (itr != m_ground->m_contact_forces.end())
        return itr->second;

    TerrainForce frc;
    frc.point = body->GetPos();
    frc.force = ChVector<>(0, 0, 0);
    frc.moment = ChVector<>(0, 0, 0);
    return frc;
}


// Set user-supplied callback for evaluating location-dependent soil parameters
void SCMDeformableTerrain::RegisterSoilParametersCallback(std::shared_ptr<SoilParametersCallback> cb) {
    m_ground->m_soil_fun = cb;
}










std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> z_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches);

bool checkRepeatVertices(ChVector<> v_temp, std::vector<ChVector<>> arr_temp);

bool CheckIdxRepeat(int target,std::vector<int> cut_cross);


int SearchVertexIdx(ChVector<> taget_vertex , std::vector<ChVector<>> vertices);

std::vector<int> GetVertexNeighbour(ChVector<> vertex, std::vector<ChSubGridMeshConnected> subMesh, std::vector<ChVector<>> vertices);

SCMDeformableSoilGrid::SCMDeformableSoilGrid(ChSystem* system, bool visualization_mesh) {
    this->SetSystem(system);

    if (visualization_mesh) {
        // Need a way to split rect mesh to tri mesh
        // Create the default mesh asset
        // m_color = std::shared_ptr<ChColorAsset>(new ChColorAsset);
        // m_color->SetColor(ChColor(0.3f, 0.3f, 0.3f));
        // this->AddAsset(m_color);

        // this->AddAsset(m_trimesh_shape);
        // m_trimesh_shape->SetWireframe(true);
    }

    do_bulldozing = false;
    bulldozing_flow_factor = 1.2;
    bulldozing_erosion_angle = 40;
    bulldozing_erosion_n_iterations = 3;
    bulldozing_erosion_n_propagations = 10;

    do_refinement = false;
    refinement_level = 1;

    // Default soil parameters
    m_Bekker_Kphi = 2e6;
    m_Bekker_Kc = 0;
    m_Bekker_n = 1.1;
    m_Mohr_cohesion = 50;
    m_Mohr_friction = 20;
    m_Janosi_shear = 0.01;
    m_elastic_K = 50000000;
    m_damping_R = 0;

   // plot_type = SCMDeformableTerrain::PLOT_NONE;
    plot_v_min = 0;
    plot_v_max = 0.2;

    test_high_offset = 0.1;
    test_low_offset = 0.5;

    last_t = 0;

    m_moving_patch = false;
}



void SCMDeformableSoilGrid::Initialize(std::shared_ptr<ChGridMeshConnected> Grid){
    m_grid_shape = Grid;
}






void SCMDeformableSoilGrid::ComputeInternalForces(){
    m_timer_calc_areas.reset();
    m_timer_ray_casting.reset();
    m_timer_refinement.reset();
    m_timer_bulldozing.reset();
    m_timer_visualization.reset();

    this->GetLoadList().clear();
    m_contact_forces.clear();

    // If enabled, update the extent of the moving patches (no ray-hit tests performed outside)
    if (m_moving_patch) {
        for (auto& p : m_patches) {
            // GET THE LOC OF THE BODY???????
            ChVector<> center_abs = p.m_body->GetFrame_REF_to_abs().TransformPointLocalToParent(p.m_point);
            ChVector<> center_loc = plane.TransformPointParentToLocal(center_abs);

            p.m_min.x() = center_loc.x() - p.m_dim.x() / 2;
            p.m_min.y() = center_loc.y() - p.m_dim.y() / 2;
            p.m_max.x() = center_loc.x() + p.m_dim.x() / 2;
            p.m_max.y() = center_loc.y() + p.m_dim.y() / 2;
        }
    }

    // Update Active SubMesh indexes
    std::vector<ChSubGridMeshConnected> subMesh= m_grid_shape->getSubGridData();

    std::vector<double> x_cut = m_grid_shape->x_cut_Arr;
    std::vector<double> z_cut = m_grid_shape->z_cut_Arr;

    std::vector<int> activeSubMesh_idx = FindActiveSubMeshIdx(x_cut, z_cut, subMesh, m_patches);


    // Get all vertices from active sub-mesh
    std::vector<ChVector<>> vertices;
    for(int i = 0; i < activeSubMesh_idx.size(); i++){
        std::vector<ChVector<>> sub_vertices = subMesh[activeSubMesh_idx[i]].getAllVertices(); 
        for(int j = 0; j<sub_vertices.size() ; j++){
            if(checkRepeatVertices(sub_vertices[j],vertices) == false){
                vertices.push_back(sub_vertices[j]);
            }
        }
    }


    // Update Connection information
    std::vector<std::vector<int>> vertices_connection;
    for(int i = 0; i < vertices.size(); i++){
        std::vector<int> neighbour_buffer = GetVertexNeighbour(vertices[i], subMesh, vertices);
        vertices_connection[i] = neighbour_buffer;
    }

    // Confirm how to structure this ?????
    // If from mesh use normal provided by mesh ???????
    ChVector<> N = plane.TransformDirectionLocalToParent(ChVector<>(0,1,0));

    // Loop through all vertices.
    // - set default SCM quantities (in case no ray-hit)
    // - skip vertices outside moving patch (if option enabled)
    // - cast ray and record result in a map (key: vertex index)
    // - initialize patch id to -1 (not set)
    struct HitRecord {
        ChContactable* contactable;  // pointer to hit object
        ChVector<> abs_point;        // hit point, expressed in global frame
        int patch_id;                // index of associated patch id
    };
    std::unordered_map<int, HitRecord> hits;

    for (int i = 0; i < vertices.size(); ++i) {
        auto v = plane.TransformParentToLocal(vertices[i]);

        // Initialize SCM quantities at current vertex
        p_sigma[i] = 0;
        p_sinkage_elastic[i] = 0;
        p_step_plastic_flow[i] = 0;
        p_erosion[i] = false;
        p_level[i] = v.z();
        p_hit_level[i] = 1e9;

        // Skip vertices outside any moving patch
        if (m_moving_patch) {
            bool outside = true;
            for (auto& p : m_patches) {
                if (v.x() >= p.m_min.x() && v.x() <= p.m_max.x() && v.y() >= p.m_min.y() && v.y() <= p.m_max.y()) {
                    outside = false;  // vertex in current patch
                    break;            // stop checking further patches
                }
            }
            if (outside)  // vertex outside all patches
                continue;
        }

        // Perform ray casting from current vertex
        collision::ChCollisionSystem::ChRayhitResult mrayhit_result;
        ChVector<> to = vertices[i] + N * test_high_offset;
        ChVector<> from = to - N * test_low_offset;
        this->GetSystem()->GetCollisionSystem()->RayHit(from, to, mrayhit_result);
        m_num_ray_casts++;
        if (mrayhit_result.hit) {
            HitRecord record = {mrayhit_result.hitModel->GetContactable(), mrayhit_result.abs_hitPoint, -1};
            hits.insert(std::make_pair(i, record));
        }


   // Loop through all hit vertices and determine to which contact patch they belong.
    // We use here the vertices_connection map (from a vertex to its adjacent vertices) which is
    // set up at initialization and updated when the mesh is refined (if refinement is enabled).
    // Use a queue-based flood-filling algorithm.
    int num_patches = 0;
    for (auto& h : hits) {
        int i = h.first;
        if (h.second.patch_id != -1)                               // move on if vertex already assigned to a patch
            continue;                                              //
        std::queue<int> todo;                                      //
        h.second.patch_id = num_patches++;                         // assign this vertex to a new patch
        todo.push(i);                                              // add vertex to end of queue
        while (!todo.empty()) {                                    //
            auto crt = hits.find(todo.front());                    // current vertex is first element in queue
            todo.pop();                                            // remove first element of queue
            auto crt_i = crt->first;                               //
            auto crt_patch = crt->second.patch_id;                 //
            for (const auto& nbr_i : vertices_connection[crt_i]) {  // loop over all neighbors
                auto nbr = hits.find(nbr_i);                       // look for neighbor in list of hit vertices
                if (nbr == hits.end())                             // move on if neighbor is not a hit vertex
                    continue;                                      //
                if (nbr->second.patch_id != -1)                    // (COULD BE REMOVED, unless we update patch area)
                    continue;                                      //
                nbr->second.patch_id = crt_patch;                  // assign neighbor to same patch
                todo.push(nbr_i);                                  // add neighbor to end of queue
            }
        }
    }

    // Collect hit vertices assigned to each patch.
    struct PatchRecord {
        std::vector<ChVector2<>> points;  // points in patch (projected on reference plane)
        double area;                      // patch area
        double perimeter;                 // patch perimeter
        double oob;                       // approximate value of 1/b
    };
    std::vector<PatchRecord> patches(num_patches);
    for (auto& h : hits) {
        ChVector<> v = plane.TransformParentToLocal(vertices[h.first]);
        patches[h.second.patch_id].points.push_back(ChVector2<>(v.x(), v.y()));
    }

    // Calculate area and perimeter of each patch.
    // Calculate approximation to Beker term 1/b.
    for (auto& p : patches) {
        utils::ChConvexHull2D ch(p.points);
        p.area = ch.GetArea();
        p.perimeter = ch.GetPerimeter();
        if (p.area < 1e-6) {
            p.oob = 0;
        } else {
            p.oob = p.perimeter / (2 * p.area);
        }
    }

    // Initialize local values for the soil parameters
    double Bekker_Kphi = m_Bekker_Kphi;
    double Bekker_Kc = m_Bekker_Kc;
    double Bekker_n = m_Bekker_n;
    double Mohr_cohesion = m_Mohr_cohesion;
    double Mohr_friction = m_Mohr_friction;
    double Janosi_shear = m_Janosi_shear;
    double elastic_K = m_elastic_K;
    double damping_R = m_damping_R;

    // Process only hit vertices
    for (auto& h : hits) {
        int i = h.first;
        ChContactable* contactable = h.second.contactable;
        const ChVector<>& abs_point = h.second.abs_point;
        int patch_id = h.second.patch_id;

        auto loc_point = plane.TransformParentToLocal(abs_point);

        /*
        if (m_soil_fun) {
            m_soil_fun->Set(loc_point.x(), loc_point.y());

            Bekker_Kphi = m_soil_fun->m_Bekker_Kphi;
            Bekker_Kc = m_soil_fun->m_Bekker_Kc;
            Bekker_n = m_soil_fun->m_Bekker_n;
            Mohr_cohesion = m_soil_fun->m_Mohr_cohesion;
            Mohr_friction = m_soil_fun->m_Mohr_friction;
            Janosi_shear = m_soil_fun->m_Janosi_shear;
            elastic_K = m_soil_fun->m_elastic_K;
            damping_R = m_soil_fun->m_damping_R;
        }
        */

        p_hit_level[i] = loc_point.y();
        double p_hit_offset = -p_hit_level[i] + p_level_initial[i];

        p_speeds[i] = contactable->GetContactPointSpeed(vertices[i]);

        ChVector<> T = -p_speeds[i];
        T = plane.TransformDirectionParentToLocal(T);
        double Vn = -T.y();
        T.y() = 0;
        T = plane.TransformDirectionLocalToParent(T);
        T.Normalize();

        // Compute i-th force:
        ChVector<> Fn;
        ChVector<> Ft;

        // Elastic try:
        p_sigma[i] = elastic_K * (p_hit_offset - p_sinkage_plastic[i]);

        // Handle unilaterality:
        if (p_sigma[i] < 0) {
            p_sigma[i] = 0;
        } else {
            // add compressive speed-proportional damping
            ////if (Vn < 0) {
            ////    p_sigma[i] += -Vn * this->damping_R;
            ////}

            p_sinkage[i] = p_hit_offset;
            p_level[i] = p_hit_level[i];

            // Accumulate shear for Janosi-Hanamoto
            p_kshear[i] += Vdot(p_speeds[i], -T) * GetSystem()->GetStep();

            // Plastic correction:
            if (p_sigma[i] > p_sigma_yeld[i]) {
                // Bekker formula
                p_sigma[i] = (patches[patch_id].oob * Bekker_Kc + Bekker_Kphi) * pow(p_sinkage[i], Bekker_n);
                p_sigma_yeld[i] = p_sigma[i];
                double old_sinkage_plastic = p_sinkage_plastic[i];
                p_sinkage_plastic[i] = p_sinkage[i] - p_sigma[i] / elastic_K;
                p_step_plastic_flow[i] = (p_sinkage_plastic[i] - old_sinkage_plastic) / GetSystem()->GetStep();
            }

            p_sinkage_elastic[i] = p_sinkage[i] - p_sinkage_plastic[i];

            // add compressive speed-proportional damping (not clamped by pressure yield)
            ////if (Vn < 0) {
            p_sigma[i] += -Vn * damping_R;
            ////}

            // Mohr-Coulomb
            double tau_max = Mohr_cohesion + p_sigma[i] * tan(Mohr_friction * CH_C_DEG_TO_RAD);

            // Janosi-Hanamoto
            p_tau[i] = tau_max * (1.0 - exp(-(p_kshear[i] / Janosi_shear)));

            Fn = N * p_area[i] * p_sigma[i];
            Ft = T * p_area[i] * p_tau[i];

            if (ChBody* rigidbody = dynamic_cast<ChBody*>(contactable)) {
                // [](){} Trick: no deletion for this shared ptr, since 'rigidbody' was not a new ChBody()
                // object, but an already used pointer because mrayhit_result.hitModel->GetPhysicsItem()
                // cannot return it as shared_ptr, as needed by the ChLoadBodyForce:
                std::shared_ptr<ChBody> srigidbody(rigidbody, [](ChBody*) {});
                std::shared_ptr<ChLoadBodyForce> mload(
                    new ChLoadBodyForce(srigidbody, Fn + Ft, false, vertices[i], false));
                this->Add(mload);

                // Accumulate contact force for this rigid body.
                // The resultant force is assumed to be applied at the body COM.
                // All components of the generalized terrain force are expressed in the global frame.
                auto itr = m_contact_forces.find(contactable);
                if (itr == m_contact_forces.end()) {
                    // Create new entry and initialize generalized force.
                    ChVector<> force = Fn + Ft;
                    TerrainForce frc;
                    frc.point = srigidbody->GetPos();
                    frc.force = force;
                    frc.moment = Vcross(Vsub(vertices[i], srigidbody->GetPos()), force);
                    m_contact_forces.insert(std::make_pair(contactable, frc));
                } else {
                    // Update generalized force.
                    ChVector<> force = Fn + Ft;
                    itr->second.force += force;
                    itr->second.moment += Vcross(Vsub(vertices[i], srigidbody->GetPos()), force);
                }
            } else if (ChLoadableUV* surf = dynamic_cast<ChLoadableUV*>(contactable)) {
                // [](){} Trick: no deletion for this shared ptr
                std::shared_ptr<ChLoadableUV> ssurf(surf, [](ChLoadableUV*) {});
                std::shared_ptr<ChLoad<ChLoaderForceOnSurface>> mload(new ChLoad<ChLoaderForceOnSurface>(ssurf));
                mload->loader.SetForce(Fn + Ft);
                mload->loader.SetApplication(0.5, 0.5);  //***TODO*** set UV, now just in middle
                this->Add(mload);

                // Accumulate contact forces for this surface.
                //// TODO
            }

            // Update mesh representation
            vertices[i] = p_vertices_initial[i] - N * p_sinkage[i];

        }  // end positive contact force

    }  // end loop on ray hits

    m_timer_ray_casting.stop();


    // ----------------------------------------------
    // REFINEMENT PENDING IMPLEMENTATION
    // ----------------------------------------------


    if (do_bulldozing) {
        std::set<int> touched_vertexes;
        for (int iv = 0; iv < vertices.size(); ++iv) {
            p_id_island[iv] = 0;
            if (p_sigma[iv] > 0)
                touched_vertexes.insert(iv);
    }

    std::set<int> domain_boundaries;

    // Compute contact islands (and their displaced material) by flood-filling the mesh
    int id_island = 0;
    for (auto fillseed = touched_vertexes.begin(); fillseed != touched_vertexes.end();
        fillseed = touched_vertexes.begin()) {
        // new island:
        ++id_island;
        std::set<int> fill_front;

        std::set<int> boundary;
        int n_vert_boundary = 0;
        double tot_area_boundary = 0;

        int n_vert_island = 1;
        double tot_step_flow_island =
            p_area[*fillseed] * p_step_plastic_flow[*fillseed] * this->GetSystem()->GetStep();
        double tot_Nforce_island = p_area[*fillseed] * p_sigma[*fillseed];
        double tot_area_island = p_area[*fillseed];
        fill_front.insert(*fillseed);
        p_id_island[*fillseed] = id_island;
        touched_vertexes.erase(fillseed);
        while (fill_front.size() > 0) {
            // fill next front
            std::set<int> fill_front_2;
            for (const auto& ifront : fill_front) {
                for (const auto& ivconnect : vertices_connection[ifront]) {
                    if ((p_sigma[ivconnect] > 0) && (p_id_island[ivconnect] == 0)) {
                        ++n_vert_island;
                        tot_step_flow_island +=
                            p_area[ivconnect] * p_step_plastic_flow[ivconnect] * this->GetSystem()->GetStep();
                        tot_Nforce_island += p_area[ivconnect] * p_sigma[ivconnect];
                        tot_area_island += p_area[ivconnect];
                        fill_front_2.insert(ivconnect);
                        p_id_island[ivconnect] = id_island;
                        touched_vertexes.erase(ivconnect);
                    } else if ((p_sigma[ivconnect] == 0) && (p_id_island[ivconnect] <= 0) &&
                                (p_id_island[ivconnect] != -id_island)) {
                        ++n_vert_boundary;
                        tot_area_boundary += p_area[ivconnect];
                        p_id_island[ivconnect] = -id_island;  // negative to mark as boundary
                        boundary.insert(ivconnect);
                    }
                }
            }
                // advance to next front
            fill_front = fill_front_2;
        }
            ////GetLog() << " island " << id_island << " flow volume =" << tot_step_flow_island << " N force=" <<
            /// tot_Nforce_island << "\n";

            // Raise the boundary because of material flow (it gives a sharp spike around the
            // island boundary, but later we'll use the erosion algorithm to smooth it out)

        for (const auto& ibv : boundary) {
            double d_y = bulldozing_flow_factor *
                         ((p_area[ibv] / tot_area_boundary) * (1 / p_area[ibv]) * tot_step_flow_island);
            double clamped_d_y = d_y;  // ChMin(d_y, ChMin(p_hit_level[ibv]-p_level[ibv], test_high_offset) );
            if (d_y > p_hit_level[ibv] - p_level[ibv]) {
                p_massremainder[ibv] += d_y - (p_hit_level[ibv] - p_level[ibv]);
                clamped_d_y = p_hit_level[ibv] - p_level[ibv];
            }
            p_level[ibv] += clamped_d_y;
            p_level_initial[ibv] += clamped_d_y;
            vertices[ibv] += N * clamped_d_y;
            p_vertices_initial[ibv] += N * clamped_d_y;
        }

        domain_boundaries.insert(boundary.begin(), boundary.end());

    }  // end for islands

    // Erosion domain area select, by topologically dilation of all the
    // boundaries of the islands:
    std::set<int> domain_erosion = domain_boundaries;
    for (const auto& ie : domain_boundaries)
        p_erosion[ie] = true;
    std::set<int> front_erosion = domain_boundaries;
    for (int iloop = 0; iloop < 10; ++iloop) {
        std::set<int> front_erosion2;
        for (const auto& is : front_erosion) {
            for (const auto& ivconnect : vertices_connection[is]) {
                if ((p_id_island[ivconnect] == 0) && (p_erosion[ivconnect] == 0)) {
                    front_erosion2.insert(ivconnect);
                    p_erosion[ivconnect] = true;
                }
            }
        }
        domain_erosion.insert(front_erosion2.begin(), front_erosion2.end());
        front_erosion = front_erosion2;
    }
    // Erosion smoothing algorithm on domain
    for (int ismo = 0; ismo < 3; ++ismo) {
        for (const auto& is : domain_erosion) {
            for (const auto& ivc : vertices_connection[is]) {
                ChVector<> vis = this->plane.TransformParentToLocal(vertices[is]);
                // flow remainder material
                if (true) {
                    if (p_massremainder[is] > p_massremainder[ivc]) {
                        double clamped_d_y_i;
                        double clamped_d_y_c;

                        // if i higher than c: clamp c upward correction as it might invalidate
                        // the ceiling constraint, if collision is nearby
                        double d_y_c = (p_massremainder[is] - p_massremainder[ivc]) *
                                        (1 / (double)vertices_connection[is].size()) * p_area[is] /
                                        (p_area[is] + p_area[ivc]);
                        clamped_d_y_c = d_y_c;
                        if (d_y_c > p_hit_level[ivc] - p_level[ivc]) {
                            p_massremainder[ivc] += d_y_c - (p_hit_level[ivc] - p_level[ivc]);
                            clamped_d_y_c = p_hit_level[ivc] - p_level[ivc];
                        }
                        double d_y_i = -d_y_c * p_area[ivc] / p_area[is];
                        clamped_d_y_i = d_y_i;
                        if (p_massremainder[is] > -d_y_i) {
                            p_massremainder[is] -= -d_y_i;
                            clamped_d_y_i = 0;
                        } else if ((p_massremainder[is] < -d_y_i) && (p_massremainder[is] > 0)) {
                            p_massremainder[is] = 0;
                            clamped_d_y_i = d_y_i + p_massremainder[is];
                        }

                        // correct vertexes
                        p_level[ivc] += clamped_d_y_c;
                        p_level_initial[ivc] += clamped_d_y_c;
                        vertices[ivc] += N * clamped_d_y_c;
                        p_vertices_initial[ivc] += N * clamped_d_y_c;

                        p_level[is] += clamped_d_y_i;
                        p_level_initial[is] += clamped_d_y_i;
                        vertices[is] += N * clamped_d_y_i;
                        p_vertices_initial[is] += N * clamped_d_y_i;
                    }
                }
                // smooth
                if (p_sigma[ivc] == 0) {
                    ChVector<> vic = this->plane.TransformParentToLocal(vertices[ivc]);
                    ChVector<> vdist = vic - vis;
                    vdist.z() = 0;
                    double ddist = vdist.Length();
                    double dy = p_level[is] + p_massremainder[is] - p_level[ivc] - p_massremainder[ivc];
                    double dy_lim = ddist * tan(bulldozing_erosion_angle * CH_C_DEG_TO_RAD);
                    if (fabs(dy) > dy_lim) {
                        double clamped_d_y_i;
                        double clamped_d_y_c;
                        if (dy > 0) {
                            // if i higher than c: clamp c upward correction as it might invalidate
                            // the ceiling constraint, if collision is nearby
                            double d_y_c = (fabs(dy) - dy_lim) * (1 / (double)vertices_connection[is].size()) *
                                            p_area[is] / (p_area[is] + p_area[ivc]);
                            clamped_d_y_c = d_y_c;  // clamped_d_y_c = ChMin(d_y_c, p_hit_level[ivc]-p_level[ivc] );
                            if (d_y_c > p_hit_level[ivc] - p_level[ivc]) {
                                p_massremainder[ivc] += d_y_c - (p_hit_level[ivc] - p_level[ivc]);
                                clamped_d_y_c = p_hit_level[ivc] - p_level[ivc];
                            }
                                double d_y_i = -d_y_c * p_area[ivc] / p_area[is];
                                clamped_d_y_i = d_y_i;
                                if (p_massremainder[is] > -d_y_i) {
                                    p_massremainder[is] -= -d_y_i;
                                    clamped_d_y_i = 0;
                                } else if ((p_massremainder[is] < -d_y_i) && (p_massremainder[is] > 0)) {
                                    p_massremainder[is] = 0;
                                    clamped_d_y_i = d_y_i + p_massremainder[is];
                                }
                            } else {
                                // if c higher than i: clamp i upward correction as it might invalidate
                                // the ceiling constraint, if collision is nearby
                                double d_y_i = (fabs(dy) - dy_lim) * (1 / (double)vertices_connection[is].size()) *
                                               p_area[is] / (p_area[is] + p_area[ivc]);
                                clamped_d_y_i = d_y_i;
                                if (d_y_i > p_hit_level[is] - p_level[is]) {
                                    p_massremainder[is] += d_y_i - (p_hit_level[is] - p_level[is]);
                                    clamped_d_y_i = p_hit_level[is] - p_level[is];
                                }
                                double d_y_c = -d_y_i * p_area[is] / p_area[ivc];
                                clamped_d_y_c = d_y_c;
                                if (p_massremainder[ivc] > -d_y_c) {
                                    p_massremainder[ivc] -= -d_y_c;
                                    clamped_d_y_c = 0;
                                } else if ((p_massremainder[ivc] < -d_y_c) && (p_massremainder[ivc] > 0)) {
                                    p_massremainder[ivc] = 0;
                                    clamped_d_y_c = d_y_c + p_massremainder[ivc];
                                }
                            }

                            // correct vertexes
                            p_level[ivc] += clamped_d_y_c;
                            p_level_initial[ivc] += clamped_d_y_c;
                            vertices[ivc] += N * clamped_d_y_c;
                            p_vertices_initial[ivc] += N * clamped_d_y_c;

                            p_level[is] += clamped_d_y_i;
                            p_level_initial[is] += clamped_d_y_i;
                            vertices[is] += N * clamped_d_y_i;
                            p_vertices_initial[is] += N * clamped_d_y_i;
                        }
                    }
                }
            }
        }

    }  // end bulldozing flow

    m_timer_bulldozing.stop();


}


}

void SCMDeformableSoilGrid::SetupAuxData(){
    m_timer_calc_areas.reset();
    m_timer_ray_casting.reset();
    m_timer_refinement.reset();
    m_timer_bulldozing.reset();
    m_timer_visualization.reset();


}




std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> z_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches){
    std::vector<int> returnBuffer;

    for(int a = 0; a<patches.size();a++){
        std::vector<int> x_cut_cross;
        std::vector<int> z_cut_cross;

        int sub_per_side = sqrt(subMesh.size());

        for(int i = 0; i<sub_per_side;i++){
            if(x_cut[i]>patches[a].m_min.x() && x_cut[i]<patches[a].m_max.x()){
                if(CheckIdxRepeat(i,x_cut_cross)==false){
                    x_cut_cross.push_back(i);
                }
                
        }

        if(z_cut[i]>patches[a].m_min.y() && z_cut[i]<patches[a].m_max.y()){
                if(CheckIdxRepeat(i,z_cut_cross)==false){
                    z_cut_cross.push_back(i);
                }
            }
        }

        x_cut_cross.push_back(x_cut_cross[x_cut_cross.size()-1]+1);
        z_cut_cross.push_back(z_cut_cross[z_cut_cross.size()-1]+1);

        for(int i = 0; i<x_cut_cross.size();i++){
            for(int j = 0; j<z_cut_cross.size();j++){
                returnBuffer.push_back(x_cut_cross[i]*sub_per_side+z_cut_cross[j]);
            }
        }
    }
    return returnBuffer;
}




bool CheckIdxRepeat(int target,std::vector<int> cut_cross){
    for(int i = 0; i<cut_cross.size();i++){
        if(target == cut_cross[i]){
            return true;
        }
    }
    return false;
}

int SearchVertexIdx(ChVector<> target_vertex , std::vector<ChVector<>> vertices){
    for(int i = 0; i<vertices.size();i++){
        if(target_vertex == vertices[i]){
            return i;
        }
    }

    return -1;
}


bool checkRepeatVertices(ChVector<> v_temp, std::vector<ChVector<>> arr_temp){
    for (int i = 0; i<arr_temp.size(); i++){
        if(arr_temp[i].x() == v_temp.x() && arr_temp[i].y() == v_temp.y() && arr_temp[i].z() == v_temp.z()){
            return true;
        }
        
    }
    return false;
}

std::vector<int> GetVertexNeighbour(ChVector<> vertex, std::vector<ChSubGridMeshConnected> subMesh, std::vector<ChVector<>> vertices)
{
    std::vector<int> return_buffer;
    for (int i = 0 ; i<subMesh.size();i++){
        std::vector<ChGridElement> ele_buffer = subMesh[i].getEleArr();
        for(int j = 0;j<ele_buffer.size();j++){
            if(ele_buffer[j].p1 == vertex){
                int idx_p2 = SearchVertexIdx(ele_buffer[j].p2,vertices);
                int idx_p3 = SearchVertexIdx(ele_buffer[j].p3,vertices);
                int idx_p4 = SearchVertexIdx(ele_buffer[j].p4,vertices);
                if(CheckIdxRepeat(idx_p2, return_buffer)){
                    return_buffer.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p3, return_buffer)){
                    return_buffer.push_back(idx_p3);
                }
                if(CheckIdxRepeat(idx_p4, return_buffer)){
                    return_buffer.push_back(idx_p4);
                }
            }

            if(ele_buffer[j].p2 == vertex){
                int idx_p1 = SearchVertexIdx(ele_buffer[j].p1,vertices);
                int idx_p3 = SearchVertexIdx(ele_buffer[j].p3,vertices);
                int idx_p4 = SearchVertexIdx(ele_buffer[j].p4,vertices);
                if(CheckIdxRepeat(idx_p1, return_buffer)){
                    return_buffer.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p3, return_buffer)){
                    return_buffer.push_back(idx_p3);
                }
                if(CheckIdxRepeat(idx_p4, return_buffer)){
                    return_buffer.push_back(idx_p4);
                }
            }

            if(ele_buffer[j].p3 == vertex){
                int idx_p1 = SearchVertexIdx(ele_buffer[j].p1,vertices);
                int idx_p2 = SearchVertexIdx(ele_buffer[j].p2,vertices);
                int idx_p4 = SearchVertexIdx(ele_buffer[j].p4,vertices);
                if(CheckIdxRepeat(idx_p1, return_buffer)){
                    return_buffer.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p2, return_buffer)){
                    return_buffer.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p4, return_buffer)){
                    return_buffer.push_back(idx_p4);
                }
            }

            if(ele_buffer[j].p4 == vertex){
                int idx_p1 = SearchVertexIdx(ele_buffer[j].p1,vertices);
                int idx_p2 = SearchVertexIdx(ele_buffer[j].p2,vertices);
                int idx_p3 = SearchVertexIdx(ele_buffer[j].p3,vertices);
                if(CheckIdxRepeat(idx_p1, return_buffer)){
                    return_buffer.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p2, return_buffer)){
                    return_buffer.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p3, return_buffer)){
                    return_buffer.push_back(idx_p3);
                }
            }
        }
    }

    return return_buffer;
}

}  // end namespace vehicle
}  // end namespace chrono
