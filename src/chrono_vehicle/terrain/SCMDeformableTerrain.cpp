
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

bool ChGridElement::ContainVertex(ChVector<double> vertex){
    if(p1 == vertex || p2 == vertex || p3 == vertex || p4 == vertex){
        return true;
    }else{
        return false;
    }
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


void ChSubGridMeshConnected::Update(ChVector<double> org, ChVector<double> new_vec){
    int idx_buff = -1;
    for(int i = 0; i<eleArr.size();i++){
        if(org==eleArr[i].p1){
            eleArr[i].p1 = new_vec;
        }
        if(org==eleArr[i].p2){
            eleArr[i].p2 = new_vec;
        }
        if(org==eleArr[i].p3){
            eleArr[i].p3 = new_vec;
        }
        if(org==eleArr[i].p4){
            eleArr[i].p4 = new_vec;
        }
    }
}

// This function performs first order refinement only
void ChSubGridMeshConnected::Refine(ChVector<double> target_vertex){
    std::vector<ChGridElement> refine_ele;
    std::vector<int> refine_ele_idx; 
    for(int i = 0; i<eleArr.size();i++){
        if(eleArr[i].p1 == target_vertex){
            refine_ele.push_back(eleArr[i]);
            refine_ele_idx.push_back(i);
            continue;
        }

        if(eleArr[i].p2 == target_vertex){
            refine_ele.push_back(eleArr[i]);
            refine_ele_idx.push_back(i);
            continue;
        }

        if(eleArr[i].p3 == target_vertex){
            refine_ele.push_back(eleArr[i]);
            refine_ele_idx.push_back(i);
            continue;
        }

        if(eleArr[i].p4 == target_vertex){
            refine_ele.push_back(eleArr[i]);
            refine_ele_idx.push_back(i);
            continue;
        }
    }

    for(int i = 0; i<refine_ele.size();i++){
        ChVector<double> p1_buff = refine_ele[i].p1;
        ChVector<double> p2_buff = refine_ele[i].p2;
        ChVector<double> p3_buff = refine_ele[i].p3;
        ChVector<double> p4_buff = refine_ele[i].p4;
        ChVector<double> normal = refine_ele[i].n;

        ChVector<double> cent_p1_p2 = (p1_buff + p2_buff)/2;
        ChVector<double> cent_p2_p3 = (p2_buff + p3_buff)/2;
        ChVector<double> cent_p3_p4 = (p3_buff + p4_buff)/2;
        ChVector<double> cent_p4_p1 = (p4_buff + p1_buff)/2;
        ChVector<double> cent = (p1_buff + p2_buff + p3_buff + p4_buff)/4;


        ChGridElement r1_buff(p1_buff, cent_p1_p2, cent, cent_p4_p1, normal);
        ChGridElement r2_buff(cent_p1_p2, p2_buff, cent_p2_p3, cent, normal);
        ChGridElement r3_buff(cent, cent_p2_p3, p3_buff, cent_p3_p4, normal);
        ChGridElement r4_buff(cent_p4_p1, cent, cent_p3_p4, p4_buff, normal);

        eleArr.push_back(r1_buff);
        eleArr.push_back(r2_buff);
        eleArr.push_back(r3_buff);
        eleArr.push_back(r4_buff);

        eleArr.erase(eleArr.begin() + refine_ele_idx[i] - i);
    }
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
    double y_tot_dis = y_max - y_min;
    double z_tot_dis = z_max - z_min;

    std::cout<<"x_tot_dis: "<<x_tot_dis<<std::endl;
    std::cout<<"y_tot_dis: "<<y_tot_dis<<std::endl;
    std::cout<<"z_tot_dis: "<<z_tot_dis<<std::endl;

    std::vector<double> x_cut;
    std::vector<double> y_cut;

    for(int i = 0; i<sub_on_side; i++){
        x_cut.push_back(x_min+(x_tot_dis/sub_on_side)*(i+1));
        y_cut.push_back(y_min+(y_tot_dis/sub_on_side)*(i+1));
    }

    x_cut_Arr = x_cut;
    y_cut_Arr = y_cut;

    // add GridElements to SubMesh and construct a connected Grid Mesh structure
    for(int i = 0; i<x_cut.size();i++){
        for(int j = 0; j<y_cut.size();j++){
            ChSubGridMeshConnected subTemp;

            for(int k = 0;k<grid_ele.size();k++){
                if(grid_ele[k].Baricenter().x()<x_cut[i] && grid_ele[k].Baricenter().y()<y_cut[j]){
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
        for(int j = 0; j<returnSubArr.size();j++){
            if(checkRepeat(returnSubArr[j],returnArr)==false){
                returnArr.push_back(returnSubArr[j]);
            }
        }
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

const std::shared_ptr<ChTriangleMeshShape> ChGridMeshConnected::GetVisMesh(){
    std::shared_ptr<ChTriangleMeshShape> return_buff = std::shared_ptr<ChTriangleMeshShape>(new ChTriangleMeshShape());
    std::shared_ptr<geometry::ChTriangleMeshConnected> convert_trimesh 
    = std::shared_ptr<geometry::ChTriangleMeshConnected>(new geometry::ChTriangleMeshConnected());

    std::vector<ChVector<double>> ret_vertices;
    std::vector<ChVector<int>> ret_faces;
    
    for(int i = 0; i<subArr.size();i++){
        std::vector<ChGridElement> eleArr_buff = subArr[i].getEleArr();
        for(int j = 0; j<eleArr_buff.size();j++){
            ChVector<double> vertex_1 = eleArr_buff[j].p1;
            ChVector<double> vertex_2 = eleArr_buff[j].p2;
            ChVector<double> vertex_3 = eleArr_buff[j].p3;
            ChVector<double> vertex_4 = eleArr_buff[j].p4;

            int idx_1 = -1;
            int idx_2 = -1;
            int idx_3 = -1;
            int idx_4 = -1;

            if(checkRepeat(vertex_1, ret_vertices)==false){
                ret_vertices.push_back(vertex_1);
                idx_1 = ret_vertices.size()-1;
            }else{
                for(int a = 0; a<ret_vertices.size();a++){
                    if(ret_vertices[a]==vertex_1){
                        idx_1 = a;
                    }
                }
            }

            if(checkRepeat(vertex_2, ret_vertices)==false){
                ret_vertices.push_back(vertex_2);
                idx_2 = ret_vertices.size()-1;
            }else{
                for(int a = 0; a<ret_vertices.size();a++){
                    if(ret_vertices[a]==vertex_2){
                        idx_2 = a;
                    }
                }
            }


            if(checkRepeat(vertex_3, ret_vertices)==false){
                ret_vertices.push_back(vertex_3);
                idx_3 = ret_vertices.size()-1;
            }else{
                for(int a = 0; a<ret_vertices.size();a++){
                    if(ret_vertices[a]==vertex_3){
                        idx_3 = a;
                    }
                }
            }

            if(checkRepeat(vertex_4, ret_vertices)==false){
                ret_vertices.push_back(vertex_4);
                idx_4 = ret_vertices.size()-1;
            }else{
                for(int a = 0; a<ret_vertices.size();a++){
                    if(ret_vertices[a]==vertex_4){
                        idx_4 = a;
                    }
                }
            }

            ChVector<int> tri_face1(idx_1, idx_2, idx_3);
            ChVector<int> tri_face2(idx_3, idx_4, idx_1);

            ret_faces.push_back(tri_face1);
            ret_faces.push_back(tri_face2);
        }
    }


    convert_trimesh->m_vertices = ret_vertices;
    convert_trimesh->m_face_v_indices = ret_faces;


    return_buff->SetMesh(convert_trimesh);

    return return_buff;
}

void ChGridMeshConnected::Update(ChVector<double> org, ChVector<double> new_vec, int submesh_idx){
    subArr[submesh_idx].Update(org, new_vec);
}

void ChGridMeshConnected::Refine(ChVector<double> target_vertex, int submesh_idx){
    subArr[submesh_idx].Refine(target_vertex);
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
                                    std::vector<double> y_cut, 
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
void SCMDeformableTerrain::Initialize(double height, double sizeX, double sizeY, int divX, int divY, int sub_per_side)
{
    std::vector<double> x_cut;
    std::vector<double> y_cut;

    double dx = sizeX / divX;
    double dy = sizeY / divY;

    std::vector<ChGridElement> temp_store_buffer;

    for(double i = 0; i<divX; i++){
        for(double j = 0; j<divY; j++){
            ChVector<double> p1;
            p1.Set(0+i*dx, 0+j*dy, height);
            ChVector<double> p2;
            p2.Set(0+(i+1)*dx, 0+j*dy, height);
            ChVector<double> p3;
            p3.Set(0+i*dx, 0+(j+1)*dy, height);
            ChVector<double> p4;
            p4.Set(0+(i+1)*dx, 0+(j+1)*dy, height);
            ChVector<double> n;
            n.Set(0, 0, 1);
            temp_store_buffer.push_back(ChGridElement(p1, p2, p3, p4,n));
        }
    }
    //Grid = new ChGridMeshConnected();
    std::cout<<"temp_store_buffer: "<<temp_store_buffer.size()<<std::endl;
    Grid->initializeData(temp_store_buffer, sub_per_side);
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

// this function might not be needed anymore
void SCMDeformableTerrain::SetAutomaticRefinementResolution(double lv) {
    m_ground->refinement_level = lv;
}

// this function might not be needed anymore
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
    m_ground->GetMesh();
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
                                    std::vector<double> y_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches);

bool checkRepeatVertices(ChVector<> v_temp, std::vector<ChVector<>> arr_temp);

bool CheckIdxRepeat(int target,std::vector<int> cut_cross);


int SearchVertexIdx(ChVector<> taget_vertex , std::vector<ChVector<>> vertices);

std::vector<int> GetVertexNeighbour(ChVector<> vertex, std::vector<ChSubGridMeshConnected> subMesh, std::vector<ChVector<>> vertices);

SCMDeformableSoilGrid::SCMDeformableSoilGrid(ChSystem* system, bool visualization_mesh) {
    this->SetSystem(system);
    m_trimesh_shape = std::shared_ptr<ChTriangleMeshShape>(new ChTriangleMeshShape);

    if (visualization_mesh) {
        // Need a way to split rect mesh to tri mesh
        // Create the default mesh asset
         m_color = std::shared_ptr<ChColorAsset>(new ChColorAsset);
         m_color->SetColor(ChColor(0.3f, 0.3f, 0.3f));
         this->AddAsset(m_color);

         this->AddAsset(m_trimesh_shape);
         m_trimesh_shape->SetWireframe(true);
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

    std::vector<double> y_cut = m_grid_shape->y_cut_Arr;

    // indexes of submesh(es) which are active
    std::vector<int> activeSubMesh_idx = FindActiveSubMeshIdx(x_cut, y_cut, subMesh, m_patches);

    // size_cut point
    std::vector<int> activeSubMesh_size_buffer;


    // Get all vertices from active sub-mesh
    vertices.clear();
    for(int i = 0; i < activeSubMesh_idx.size(); i++){
        std::vector<ChVector<>> sub_vertices = subMesh[activeSubMesh_idx[i]].getAllVertices(); 
        for(int j = 0; j<sub_vertices.size() ; j++){
            if(checkRepeatVertices(sub_vertices[j],vertices) == false){
                vertices.push_back(sub_vertices[j]);
            }
        }
        activeSubMesh_size_buffer.push_back(vertices.size()-1); //push the index cut point of that array
    }

    SetupAuxData();

    std::cout<<"test point 1"<<std::endl;


    // Update Connection information
    std::vector<std::vector<int>> vertices_connection;
    for(int i = 0; i < vertices.size(); i++){
        std::vector<int> neighbour_buffer = GetVertexNeighbour(vertices[i], subMesh, vertices);
        vertices_connection.push_back(neighbour_buffer);
    }
    for (unsigned int iv = 0; iv < vertices.size(); ++iv) {
        p_area[iv] = 0,04; //need an api to get uniform area???????
    }


    // Confirm how to structure this ?????
    // If from mesh use normal provided by mesh ???????
    //ChVector<> N = plane.TransformDirectionLocalToParent(ChVector<>(0,0,1));
    ChVector<> N = ChVector<>(0,0,1);
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
     std::cout<<"test point 2"<<std::endl;
    std::vector<int> hit_vertices_idx;
    for (int i = 0; i < vertices.size(); ++i) {
        //auto v = plane.TransformParentToLocal(vertices[i]);
        auto v = vertices[i];

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
            hit_vertices_idx.push_back(i);
        }
    }
     std::cout<<"test point 3"<<std::endl;
    // this vector stores all original vertices
    std::vector<ChVector<double>> original_vertice_hit;
    std::cout<<"hit_size: "<<hit_vertices_idx.size()<<std::endl;
    // Now we have the indexes of all hit vertices, now we need to store all these original ChVector 
    // for comparison
    for (int i = 0; i<hit_vertices_idx.size();i++){
        original_vertice_hit.push_back(vertices[i]);
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

            }

            // Update mesh representation
            vertices[i] = p_vertices_initial[i] - N * p_sinkage[i];

        }  // end positive contact force

    }  // end loop on ray hits

    m_timer_ray_casting.stop();

    // -----------------------------------------------
    // store all updated hit vertices
    // -----------------------------------------------
    std::vector<ChVector<double>> processed_vertices_hit;
    for (int i = 0; i<hit_vertices_idx.size();i++){
        processed_vertices_hit.push_back(vertices[i]);
        std::cout<<"stoer all updated hit"<<std::endl;
    }
    // ----------------------------------------------
    // ----------------------------------------------


    // ----------------------------------------------
    // update the Grid data structure
    // ----------------------------------------------
    for (int i = 0; i<hit_vertices_idx.size();i++){
           for(int j = 0; j<activeSubMesh_idx.size(); j++){
            if(hit_vertices_idx[i]>activeSubMesh_size_buffer[j]){
                m_grid_shape->Update(original_vertice_hit[hit_vertices_idx[i]],vertices[hit_vertices_idx[i]], j+1);
                break;
            }

            if(j==activeSubMesh_idx.size()-1){
                m_grid_shape->Update(original_vertice_hit[hit_vertices_idx[i]],vertices[hit_vertices_idx[i]], 0);
            }
        }
        std::cout<<"update the Grid data structure"<<std::endl;
    }
    // ----------------------------------------------
    // ----------------------------------------------








    // ----------------------------------------------
    // Rectangular Mesh Refinement
    // ----------------------------------------------
    for(int i = 0; i<hit_vertices_idx.size();i++){
        for(int j = 0; j<activeSubMesh_idx.size(); j++){
            if(hit_vertices_idx[i]>activeSubMesh_size_buffer[j]){
                m_grid_shape->Refine(vertices[hit_vertices_idx[i]], j+1);
                break;
            }

            if(j==activeSubMesh_idx.size()-1){
                m_grid_shape->Refine(vertices[hit_vertices_idx[i]], 0);
            }
        }
        
    }
    // ----------------------------------------------
    // ----------------------------------------------




}

void SCMDeformableSoilGrid::SetupAuxData(){
    p_speeds.resize(vertices.size());
    p_step_plastic_flow.resize(vertices.size());
    p_level.resize(vertices.size());
    p_level_initial.resize(vertices.size());
    p_hit_level.resize(vertices.size());
    p_sinkage.resize(vertices.size());
    p_sinkage_plastic.resize(vertices.size());
    p_sinkage_elastic.resize(vertices.size());
    p_kshear.resize(vertices.size());
    p_area.resize(vertices.size());
    p_sigma.resize(vertices.size());
    p_sigma_yeld.resize(vertices.size());
    p_tau.resize(vertices.size());
    p_massremainder.resize(vertices.size());
    p_id_island.resize(vertices.size());
    p_erosion.resize(vertices.size());

    for (int i = 0; i < vertices.size(); ++i) {
        p_level[i] = plane.TransformParentToLocal(vertices[i]).z();
        p_level_initial[i] = p_level[i];
    }
}

void SCMDeformableSoilGrid::GetMesh(){
    m_trimesh_shape = m_grid_shape->GetVisMesh();
}



bool CheckIdxRepeat(int target,std::vector<int> cut_cross);
std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> y_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches){
    std::vector<int> returnBuffer;

    for(int i = 0; i<subMesh.size();i++){
        subMesh[i].getBoundingInfo();
        std::vector<ChGridElement> ele_arr_temp = subMesh[i].getEleArr();
        for(int j = 0; j<patches.size();j++){
            
            double patch_cen_x = (patches[j].m_min.x() + patches[j].m_max.x()) / 2;
            double patch_cen_y = (patches[j].m_min.y() + patches[j].m_max.y()) / 2;

            if(patch_cen_x > subMesh[i].xmin && patch_cen_x <subMesh[i].xmax && patch_cen_y > subMesh[i].ymin && patch_cen_y < subMesh[i].ymax)
            {
                returnBuffer.push_back(i);
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
