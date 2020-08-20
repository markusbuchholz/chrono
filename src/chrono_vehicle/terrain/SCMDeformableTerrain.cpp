
#include "chrono_vehicle/terrain/SCMDeformableTerrain.h"
#include <math.h> 
#include <queue>


namespace chrono {
namespace vehicle {


bool CheckIdxRepeat(int target,std::vector<int> cut_cross){
    for(int i = 0; i<cut_cross.size();i++){
        if(target == cut_cross[i]){
            return true;
        }
    }
    return false;
}

int searchVertexIdx(ChVector<> target_vertex , std::vector<ChVector<>> vertices){
    for(int i = 0; i<vertices.size();i++){
        if(target_vertex == vertices[i]){
            return i;
        }
    }

    return -1;
}

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


std::vector<ChVector<double>> ChSubGridMeshConnected::getAllVertices_vec(){
    return vertices_vec;
}

std::vector<std::vector<int>> ChSubGridMeshConnected::getAllNeigh_vec(){
    return neighbour_map_vec;
}

std::vector<ChVector<int>> ChSubGridMeshConnected::getAllFaces(){
    return face_vec;
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




void ChSubGridMeshConnected::Update(ChVector<double> new_vec, int idx){
    vertices_vec[idx] = new_vec;
    /*for(int i = 0; i<eleArr.size();i++){

        //std::cout<<"============"<<std::endl;
        //std::cout<<"i: "<<i<<" x: "<<org.x()<<" comp x: "<<eleArr[i].p1.x()<<std::endl;
        if(org.x()==eleArr[i].p1.x() && org.y()==eleArr[i].p1.y()){
            //std::cout<<"change p1 aha old: "<<org<<" new: "<<new_vec<<std::endl;
            //eleArr[i].p1.x() = new_vec.x();
            //eleArr[i].p1.y() = new_vec.y();
            eleArr[i].p1.z() = new_vec.z();
            //std::cout<<"change p1 old: "<<org<<" new: "<<new_vec<<std::endl;
        }
        if(org.x()==eleArr[i].p2.x() && org.y()==eleArr[i].p2.y()){
            //eleArr[i].p2.x() = new_vec.x();
            //eleArr[i].p2.y() = new_vec.y();
            eleArr[i].p2.z() = new_vec.z();
            //std::cout<<"change p2 old: "<<org<<" new: "<<new_vec<<std::endl;
        }
        if(org.x()==eleArr[i].p3.x() && org.y()==eleArr[i].p3.y()){
            //eleArr[i].p3.x() = new_vec.x();
            //eleArr[i].p3.y() = new_vec.y();
            eleArr[i].p3.z() = new_vec.z();
            //std::cout<<"change p3 old: "<<org<<" new: "<<new_vec<<std::endl;
        }
        if((org.x()==eleArr[i].p4.x() && org.y()==eleArr[i].p4.y())){
            //eleArr[i].p4.x() = new_vec.x();
            //eleArr[i].p4.y() = new_vec.y();
            eleArr[i].p4.z() = new_vec.z();
            //std::cout<<"change p4 old: "<<org<<" new: "<<new_vec<<std::endl;
        }
    }*/

    
}

void ChSubGridMeshConnected::UpdateColor(ChVector<float> new_color, int idx){
    colors_vec[idx] = new_color;
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

void ChSubGridMeshConnected::InitAllVertices(){

    vertices_vec.clear();
    face_vec.clear();

    for(int i = 0; i<eleArr.size();i++){
        std::vector<ChVector<double>> vertices_temp;
        if(checkRepeat(eleArr[i].p1, vertices_vec)==false){
            vertices_temp.push_back(eleArr[i].p1);
        }

        if(checkRepeat(eleArr[i].p2, vertices_vec)==false){
            vertices_temp.push_back(eleArr[i].p2);
        }

        if(checkRepeat(eleArr[i].p3, vertices_vec)==false){
            vertices_temp.push_back(eleArr[i].p3);
        }

        if(checkRepeat(eleArr[i].p4, vertices_vec)==false){
            vertices_temp.push_back(eleArr[i].p4);
        }

        vertices_vec.insert(vertices_vec.end(),vertices_temp.begin(), vertices_temp.end());
    }

    for(int j = 0; j<eleArr.size();j++){
        int idx_1 = searchVertexIdx(eleArr[j].p1,vertices_vec);
        int idx_2 = searchVertexIdx(eleArr[j].p2,vertices_vec);
        int idx_3 = searchVertexIdx(eleArr[j].p3,vertices_vec);
        int idx_4 = searchVertexIdx(eleArr[j].p4,vertices_vec);

        ChVector<int> tri_face1(idx_1, idx_2, idx_4);
        ChVector<int> tri_face2(idx_4, idx_3, idx_1);

        face_vec.push_back(tri_face1);
        face_vec.push_back(tri_face2);
    }
}




int SearchSubVertexIdx(ChVector<> vertex, std::vector<ChVector<>> arr){
    int return_value = -1;
    for(int i = 0;i<arr.size();i++){
        if(arr[i] == vertex){
            return_value = i;
        }
    }
    return return_value;
}

void ChSubGridMeshConnected::InitVerColor(){
    for(int i = 0; i<vertices_vec.size();i++){
        colors_vec.push_back(ChVector<>(0.3f, 0.3f, 0.3f));
    }
}

std::vector<ChVector<>> ChSubGridMeshConnected::getAllColors_vec(){
    return colors_vec;
}


void ChSubGridMeshConnected::InitVerNeighMap(){
    //vertices_vec.clear();
    neighbour_map_vec.clear();
    for(int i = 0; i<vertices_vec.size();i++){
        std::vector<int> neighbour_idx;
        for(int j = 0;j<eleArr.size();j++){
            if(eleArr[j].p1 == vertices_vec[i]){
            int idx_p2 = SearchSubVertexIdx(eleArr[j].p2,vertices_vec);
            int idx_p3 = SearchSubVertexIdx(eleArr[j].p3,vertices_vec);
            int idx_p4 = SearchSubVertexIdx(eleArr[j].p4,vertices_vec);
                if(CheckIdxRepeat(idx_p2, neighbour_idx)){
                    neighbour_idx.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p3, neighbour_idx)){
                    neighbour_idx.push_back(idx_p3);
                }
                if(CheckIdxRepeat(idx_p4, neighbour_idx)){
                    neighbour_idx.push_back(idx_p4);
                }
            }

            if(eleArr[j].p2 == vertices_vec[i]){
            int idx_p1 = SearchSubVertexIdx(eleArr[j].p1,vertices_vec);
            int idx_p3 = SearchSubVertexIdx(eleArr[j].p3,vertices_vec);
            int idx_p4 = SearchSubVertexIdx(eleArr[j].p4,vertices_vec);
                if(CheckIdxRepeat(idx_p1, neighbour_idx)){
                    neighbour_idx.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p3, neighbour_idx)){
                    neighbour_idx.push_back(idx_p3);
                }
                if(CheckIdxRepeat(idx_p4, neighbour_idx)){
                    neighbour_idx.push_back(idx_p4);
                }
            }

            if(eleArr[j].p3 == vertices_vec[i]){
            int idx_p1 = SearchSubVertexIdx(eleArr[j].p1,vertices_vec);
            int idx_p2 = SearchSubVertexIdx(eleArr[j].p2,vertices_vec);
            int idx_p4 = SearchSubVertexIdx(eleArr[j].p4,vertices_vec);
                if(CheckIdxRepeat(idx_p1, neighbour_idx)){
                    neighbour_idx.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p2, neighbour_idx)){
                    neighbour_idx.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p4, neighbour_idx)){
                    neighbour_idx.push_back(idx_p4);
                }
            }

            if(eleArr[j].p4 == vertices_vec[i]){
            int idx_p1 = SearchSubVertexIdx(eleArr[j].p1,vertices_vec);
            int idx_p2 = SearchSubVertexIdx(eleArr[j].p2,vertices_vec);
            int idx_p3 = SearchSubVertexIdx(eleArr[j].p3,vertices_vec);
                if(CheckIdxRepeat(idx_p1, neighbour_idx)){
                    neighbour_idx.push_back(idx_p1);
                }
                if(CheckIdxRepeat(idx_p2, neighbour_idx)){
                    neighbour_idx.push_back(idx_p2);
                }
                if(CheckIdxRepeat(idx_p3, neighbour_idx)){
                    neighbour_idx.push_back(idx_p3);
                }
            }
        }
        neighbour_map_vec.push_back(neighbour_idx);
    }
}



std::vector<ChVector<>> ChSubGridMeshConnected::returnMeshVert(){
    return mesh_vertices;
}
std::vector<ChVector<int>> ChSubGridMeshConnected::returnMeshFace(){
    return mesh_face;
}



bool checkRepeatPtr(ChVector<> vex, std::vector<ChVector<>>&);
void ChSubGridMeshConnected::GetSubVisMesh(ChCoordsys<> plane){
    mesh_vertices.clear();
    mesh_face.clear();
    //std::cout<<"test 1dsdsd"<<std::endl;

    for(int j = 0; j<eleArr.size();j++){
        ChVector<double> vertex_1 = plane.TransformPointLocalToParent(eleArr[j].p1);
        ChVector<double> vertex_2 = plane.TransformPointLocalToParent(eleArr[j].p2);
        ChVector<double> vertex_3 = plane.TransformPointLocalToParent(eleArr[j].p3);
        ChVector<double> vertex_4 = plane.TransformPointLocalToParent(eleArr[j].p4);

        //std::cout<<"y: "<<vertex_1.y()<<std::endl;
        //std::cout<<"y: "<<vertex_2.y()<<std::endl;
        //std::cout<<"y: "<<vertex_3.y()<<std::endl;
        //std::cout<<"y: "<<vertex_4.y()<<std::endl;

        int idx_1 = -1;
        int idx_2 = -1;
        int idx_3 = -1;
        int idx_4 = -1;

        if(checkRepeatPtr(vertex_1, mesh_vertices)==false){
            mesh_vertices.push_back(vertex_1);
            idx_1 = mesh_vertices.size()-1;
        }else{
            for(int a = 0; a<mesh_vertices.size();a++){
                if(mesh_vertices[a]==vertex_1){
                    idx_1 = a;
                }
            }
        }

        if(checkRepeatPtr(vertex_2, mesh_vertices)==false){
            mesh_vertices.push_back(vertex_2);
            idx_2 = mesh_vertices.size()-1;
        }else{
            for(int a = 0; a<mesh_vertices.size();a++){
                if(mesh_vertices[a]==vertex_2){
                    idx_2 = a;
                }
            }
        }


        if(checkRepeatPtr(vertex_3, mesh_vertices)==false){
            mesh_vertices.push_back(vertex_3);
            idx_3 = mesh_vertices.size()-1;
        }else{
            for(int a = 0; a<mesh_vertices.size();a++){
                if(mesh_vertices[a]==vertex_3){
                    idx_3 = a;
                }
            }
        }

        if(checkRepeatPtr(vertex_4, mesh_vertices)==false){
            mesh_vertices.push_back(vertex_4);
            idx_4 = mesh_vertices.size()-1;
        }else{
            for(int a = 0; a<mesh_vertices.size();a++){
                if(mesh_vertices[a]==vertex_4){
                    idx_4 = a;
                }
            }
        }



        ChVector<int> tri_face1(idx_1, idx_2, idx_4);
        ChVector<int> tri_face2(idx_4, idx_3, idx_1);


        mesh_face.push_back(tri_face1);
        mesh_face.push_back(tri_face2);

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

    int total_sub = sub_on_side * sub_on_side;
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
                if(grid_ele[k].Baricenter().x()<=x_cut[i] && grid_ele[k].Baricenter().y()<=y_cut[j]){
                    subTemp.addGridElement(grid_ele[k]);
                    grid_ele.erase(grid_ele.begin()+k);
                    k=k-1;
                }
               
            }
            std::cout<<"Preprocessing SubMesh: "<<subArr.size()<<" total Submesh: "<<total_sub<<std::endl;
            addSubGridData(subTemp);
        }
    }

    for(int i = 0; i<subArr.size();i++){
        subArr[i].getBoundingInfo();
    }


    InitSubAllVertices();

    InitSubAllColors();

    InitializeSubNeighMap();
    
}

void ChGridMeshConnected::addSubGridData(ChSubGridMeshConnected subMesh){
    subArr.push_back(subMesh);
}




std::vector<ChSubGridMeshConnected> ChGridMeshConnected::getSubGridData(){
    // return a vector of sub mesh
    return subArr;
}


void ChGridMeshConnected::InitializeSubNeighMap(){
    for(int i = 0; i<subArr.size();i++){
        std::cout<<"setting up neighbour map for i: "<<i<<std::endl;
        subArr[i].InitVerNeighMap();
    }

}

void ChGridMeshConnected::InitSubAllColors(){
    for(int i = 0; i<subArr.size();i++){
        std::cout<<"setting up color vis for i: "<<i<<std::endl;
        subArr[i].InitVerColor();
    }
}

void ChGridMeshConnected::InitSubAllVertices(){
    for(int i = 0; i<subArr.size();i++){
        std::cout<<"setting up vertices for i: "<<i<<std::endl;
        subArr[i].InitAllVertices();
    }
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




void ChGridMeshConnected::GetVisMesh(std::shared_ptr<ChTriangleMeshShape> trimesh, ChCoordsys<> plane
,std::vector<int> active_sub_mesh){
    std::vector<ChVector<int>>& idx_vertices = trimesh->GetMesh()->getIndicesVertexes();
    std::vector<ChVector<>>& vertices = trimesh->GetMesh()->getCoordsVertices();
    std::vector<ChVector<float>>& colors = trimesh->GetMesh()->getCoordsColors();

    idx_vertices.clear();
    vertices.clear();
    colors.clear();
    
    int idx_size_ind = 0;

    
    for(int i = 0; i<subArr.size();i++){
        std::vector<ChVector<>> ver_buff = subArr[i].getAllVertices_vec();
        std::vector<ChVector<>> color_buff = subArr[i].getAllColors_vec();
        std::vector<ChVector<int>> face_buff = subArr[i].getAllFaces();

        for(int a = 0; a<ver_buff.size();a++){
            ver_buff[a] = plane.TransformPointLocalToParent(ver_buff[a]);
        }

        vertices.insert(vertices.end(),ver_buff.begin(),ver_buff.end());

        for(int a = 0; a<face_buff.size();a++){
            face_buff[a].x() = face_buff[a].x() + idx_size_ind;
            face_buff[a].y() = face_buff[a].y() + idx_size_ind;
            face_buff[a].z() = face_buff[a].z() + idx_size_ind;
        }

        idx_vertices.insert(idx_vertices.end(),face_buff.begin(),face_buff.end());

        colors.insert(colors.end(),color_buff.begin(),color_buff.end());

        idx_size_ind = vertices.size();

    }
    

}

void ChGridMeshConnected::InitializeMeshVis(ChCoordsys<> plane){
    for(int i = 0 ; i<subArr.size();i++){
        subArr[i].GetSubVisMesh(plane);
    }
}





bool checkRepeatPtr(ChVector<> vex, std::vector<ChVector<>>& arr){
    for(int i = 0; i<arr.size();i++){
        if(arr[i] == vex){
            return true;
        }
    }

    return false;

}

void ChGridMeshConnected::Update(ChVector<double> new_vec,int idx ,int submesh_idx){
    //if(submesh_idx==2){
        subArr[submesh_idx].Update(new_vec, idx);
        //std::cout<<"old: "<<org<<" new: "<<new_vec<<std::endl;
    //}
}

void ChGridMeshConnected::UpdateColor(ChVector<float> new_color,int idx ,int submesh_idx){
    subArr[submesh_idx].UpdateColor(new_color, idx);
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



SCMDeformableTerrain::SCMDeformableTerrain(ChSystem* system, bool visualization_mesh ) {
    m_ground = chrono_types::make_shared<SCMDeformableSoilGrid>(system, visualization_mesh);
    system->Add(m_ground);
    Grid = std::make_shared<ChGridMeshConnected>();
}

// Initialize a grid mesh by input values
void SCMDeformableTerrain::Initialize(double height, double sizeX, double sizeY, int divX, int divY, int sub_per_side, ChCoordsys<> plane)
{
    SetPlane(plane);

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

    Grid->initializeData(temp_store_buffer, sub_per_side);

    

    Grid->InitializeMeshVis(plane);

    m_ground->Initialize(Grid);

    m_ground->m_height = height;

    m_ground->area_x = sizeX/divX;
    m_ground->area_y = sizeY/divY;
}


void SCMDeformableTerrain::Initialize(const std::string& mesh_file, int sub_per_side,ChCoordsys<> plane)
{
    GridMeshLoader test;
    std::vector<ChGridElement> result = test.loadObj(mesh_file);
    Grid->initializeData(result, sub_per_side);

    SetPlane(plane);
    
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
    os << "   Vertices Setup:           " << m_ground->m_timer_vertsetup() << std::endl;
    os << "   Visualization:           " << m_ground->m_timer_visualization() << std::endl;
    os << "   Update Time:           " << m_ground->m_timer_update() << std::endl;
    os << "   Total Time:           " << m_ground->m_timer_total() << std::endl;

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

void SCMDeformableTerrain::SetColor(ChColor color) {
    if (m_ground->m_color)
        m_ground->m_color->SetColor(color);
}



std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> y_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches);

bool checkRepeatVertices(ChVector<> v_temp, std::vector<ChVector<>> arr_temp);



int SearchVertexIdx(ChVector<> taget_vertex , std::vector<ChVector<>> vertices);

std::vector<int> GetVertexNeighbour(ChVector<> vertex, std::vector<ChSubGridMeshConnected> subMesh, std::vector<ChVector<>> vertices);

SCMDeformableSoilGrid::SCMDeformableSoilGrid(ChSystem* system, bool visualization_mesh) {
    this->SetSystem(system);
    m_trimesh_shape = std::shared_ptr<ChTriangleMeshShape>(new ChTriangleMeshShape);

    if (visualization_mesh) {
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
    m_timer_update.reset();
    m_timer_refinement.reset();
    m_timer_bulldozing.reset();
    m_timer_visualization.reset();
    m_timer_total.reset();
    m_timer_vertsetup.reset();




    m_timer_total.start();

    

    this->GetLoadList().clear();
    m_contact_forces.clear();

    // If enabled, update the extent of the moving patches (no ray-hit tests performed outside)
    if (m_moving_patch) {
        for (auto& p : m_patches) {
            ChVector<> center_abs = p.m_body->GetFrame_REF_to_abs().TransformPointLocalToParent(p.m_point);
            ChVector<> center_loc = plane.TransformPointParentToLocal(center_abs);
            p.m_min.x() = center_loc.x() - p.m_dim.x() / 2;
            p.m_min.y() = center_loc.y() - p.m_dim.y() / 2;
            p.m_max.x() = center_loc.x() + p.m_dim.x() / 2;
            p.m_max.y() = center_loc.y() + p.m_dim.y() / 2;
        }
    }else {
        UpdateFixedPatch();
    }

   

    // Update Active SubMesh indexes
    std::vector<ChSubGridMeshConnected> subMesh= m_grid_shape->getSubGridData();

    std::vector<double> x_cut = m_grid_shape->x_cut_Arr;

    std::vector<double> y_cut = m_grid_shape->y_cut_Arr;

    // indexes of submesh(es) which are active
    active_sub_mesh.clear();
    m_timer_vertsetup.start();
    //std::cout<<"active_sub_mesh_length: "<<active_sub_mesh.size()<<std::endl;
    active_sub_mesh = FindActiveSubMeshIdx(x_cut, y_cut, subMesh, m_patches);
    m_timer_vertsetup.stop();
    //std::cout<<"active_sub_mesh_length after: "<<active_sub_mesh.size()<<std::endl;
    
    std::vector<int> activeSubMesh_size_buffer;

   // Get all vertices from active sub-mesh
    vertices.clear();
    //std::cout<<"vertices size before add: "<<vertices.size()<<std::endl;
    int neigh_size_ind = 0;
    std::vector<std::vector<int>> vertices_connection;
    for(int i = 0; i < active_sub_mesh.size(); i++){
        std::vector<ChVector<>> sub_vertices = subMesh[active_sub_mesh[i]].getAllVertices_vec(); 
        std::cout<<"partition: "<<active_sub_mesh[i]<<" size: "<<sub_vertices.size()<<std::endl;
        vertices.insert(vertices.end(),sub_vertices.begin(),sub_vertices.end());

        activeSubMesh_size_buffer.push_back(vertices.size()-1); //push the index cut point of that array
    
        std::vector<std::vector<int>> neighbour_buffer = subMesh[active_sub_mesh[i]].getAllNeigh_vec();
        for(int a = 0; a<neighbour_buffer.size();a++){
            for(int b = 0; b<neighbour_buffer[a].size();b++){
                neighbour_buffer[a][b] = neighbour_buffer[a][b]+neigh_size_ind;
            }
        }
        vertices_connection.insert(vertices_connection.end(),neighbour_buffer.begin(),neighbour_buffer.end());
        neigh_size_ind = vertices.size();
    }

    std::cout<<"vertices size: "<<vertices.size()<<std::endl;

    


    m_timer_calc_areas.start();
    SetupAuxData();

    for (int i = 0; i < vertices.size(); i++) {
        p_level_initial[i] = (vertices[i]).z();
        p_area[i] = area_x*area_y*10*(m_height - p_level_initial[i]+1)*(m_height - p_level_initial[i]+1); //need an api to get uniform area???????
    }

    m_timer_calc_areas.stop();

    // Confirm how to structure this ?????
    // If from mesh use normal provided by mesh ???????
    ChVector<> N = plane.TransformDirectionLocalToParent(ChVector<>(0,0,1));
    //ChVector<> N = ChVector<>(0,0,1);
    // Loop through all vertices.
    // - set default SCM quantities (in case no ray-hit)
    // - skip vertices outside moving patch (if option enabled)
    // - cast ray and record result in a map (key: vertex index)
    // - initialize patch id to -1 (not set)


    m_timer_ray_casting.start();

    struct HitRecord {
        ChContactable* contactable;  // pointer to hit object
        ChVector<> abs_point;        // hit point, expressed in global frame
        int patch_id;                // index of associated patch id
    };
    std::unordered_map<int, HitRecord> hits;
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
        
        bool outside = true;
        for (auto& p : m_patches) {

            if (v.x() >= p.m_min.x() && v.x() <= p.m_max.x() && v.y() >= p.m_min.y() && v.y() <= p.m_max.y()) {
                outside = false;  // vertex in current patch
                break;            // stop checking further patches
            }
        }
        if (outside)  // vertex outside all patches
            continue;
    

        // Perform ray casting from current vertex
        collision::ChCollisionSystem::ChRayhitResult mrayhit_result;
        //ChVector<> to = vertices[i] + N * test_high_offset;
        ChVector<> to = plane.TransformPointLocalToParent(vertices[i]+ChVector<>(0,0,1)*test_high_offset);
        ChVector<> from = to - N * test_low_offset;

        this->GetSystem()->GetCollisionSystem()->RayHit(from, to, mrayhit_result);
        m_num_ray_casts++;

        if (mrayhit_result.hit) {
            HitRecord record = {mrayhit_result.hitModel->GetContactable(), mrayhit_result.abs_hitPoint, -1};
            hits.insert(std::make_pair(i, record));
            hit_vertices_idx.push_back(i);
        }
    }


    // this vector stores all original vertices
    std::vector<ChVector<double>> original_vertice_hit;
    //std::cout<<"hit_size: "<<hit_vertices_idx.size()<<std::endl;
    // Now we have the indexes of all hit vertices, now we need to store all these original ChVector 
    // for comparison
    for (int i = 0; i<hit_vertices_idx.size();i++){
        original_vertice_hit.push_back(vertices[hit_vertices_idx[i]]);
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
            auto crt_patch = crt->second.patch_id;              //
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
    //std::cout<<"test point 7"<<std::endl;
    // Collect hit vertices assigned to each patch.
    struct PatchRecord {
        std::vector<ChVector2<>> points;  // points in patch (projected on reference plane)
        double area;                      // patch area
        double perimeter;                 // patch perimeter
        double oob;                       // approximate value of 1/b
    };
    std::vector<PatchRecord> patches(num_patches);
    for (auto& h : hits) {
        //ChVector<> v = plane.TransformParentToLocal(vertices[h.first]);
        ChVector<> v = vertices[h.first];

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
    //std::cout<<"test point 8"<<std::endl;

    // Initialize local values for the soil parameters
    double Bekker_Kphi = m_Bekker_Kphi;
    double Bekker_Kc = m_Bekker_Kc;
    double Bekker_n = m_Bekker_n;
    double Mohr_cohesion = m_Mohr_cohesion;
    double Mohr_friction = m_Mohr_friction;
    double Janosi_shear = m_Janosi_shear;
    double elastic_K = m_elastic_K;
    double damping_R = m_damping_R;
    //std::cout<<"test point 9"<<std::endl;
    // Process only hit vertices
    for (auto& h : hits) {
        int i = h.first;
        ChContactable* contactable = h.second.contactable;
        const ChVector<>& abs_point = h.second.abs_point;
        int patch_id = h.second.patch_id;
        //std::cout<<"test point 10"<<std::endl;

        auto loc_point = plane.TransformParentToLocal(abs_point);
        //auto loc_point = abs_point;
        //std::cout<<"test point 11"<<std::endl;
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

        p_hit_level[i] = loc_point.z();
        double p_hit_offset = -p_hit_level[i] + p_level_initial[i];
        //std::cout<<"p_hit_offset: "<<p_hit_offset<<std::endl;

        p_speeds[i] = contactable->GetContactPointSpeed(plane.TransformLocalToParent(vertices[i]));
        //std::cout<<"test poitn 7"<<std::endl;
        //std::cout<<"p_speeds[]: "<<p_speeds[i]<<std::endl;
        ChVector<> T = -p_speeds[i];
        //T = plane.TransformDirectionParentToLocal(T);
        //std::cout<<"T: "<<T<<std::endl;
        double Vn = -T.y();
        T.y() = 0;
        //T = plane.TransformDirectionLocalToParent(T);
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
                    //std::cout<<"force: "<<force<<std::endl;
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

            vertices[i] = p_vertices_initial[i] - ChVector<>(0,0,1) * p_sinkage[i];

        }  // end positive contact force

    }  // end loop on ray hits
    m_timer_ray_casting.stop();



    // -----------------------------------------------
    // store all updated hit vertices
    // -----------------------------------------------
    //std::cout<<"processed size: "<<hit_vertices_idx.size()<<std::endl;

    std::vector<ChVector<double>> processed_vertices_hit;

    for (int i = 0; i<hit_vertices_idx.size();i++){

        processed_vertices_hit.push_back(vertices[hit_vertices_idx[i]]);
    }


    // ----------------------------------------------
    // update the Grid data structure
    // ----------------------------------------------






    



    
    m_timer_visualization.start();
    // Create the default mesh asset
    m_color = std::shared_ptr<ChColorAsset>(new ChColorAsset);
    m_color->SetColor(ChColor(0.3f, 0.3f, 0.3f));
    this->AddAsset(m_color);
    GetMesh();
    this->AddAsset(m_trimesh_shape);
    m_trimesh_shape->SetWireframe(true);


   
    // ----------------------------------------------
    // ----------------------------------------------

    //
    // Update the visualization colors
    //

    std::vector<ChVector<float>> colors; 

    if (plot_type != SCMDeformableTerrain::PLOT_NONE) {
        colors.resize(hit_vertices_idx.size());
        for (size_t iv = 0; iv < hit_vertices_idx.size(); ++iv) {
            ChColor mcolor;
            switch (plot_type) {
                case SCMDeformableTerrain::PLOT_LEVEL:
                    mcolor = ChColor::ComputeFalseColor(p_level[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_LEVEL_INITIAL:
                    mcolor = ChColor::ComputeFalseColor(p_level_initial[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_SINKAGE:
                    mcolor = ChColor::ComputeFalseColor(p_sinkage[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_SINKAGE_ELASTIC:
                    mcolor = ChColor::ComputeFalseColor(p_sinkage_elastic[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_SINKAGE_PLASTIC:
                    mcolor = ChColor::ComputeFalseColor(p_sinkage_plastic[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_STEP_PLASTIC_FLOW:
                    mcolor = ChColor::ComputeFalseColor(p_step_plastic_flow[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_K_JANOSI:
                    mcolor = ChColor::ComputeFalseColor(p_kshear[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_PRESSURE:
                    mcolor = ChColor::ComputeFalseColor(p_sigma[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_PRESSURE_YELD:
                    mcolor = ChColor::ComputeFalseColor(p_sigma_yeld[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_SHEAR:
                    mcolor = ChColor::ComputeFalseColor(p_tau[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_MASSREMAINDER:
                    mcolor = ChColor::ComputeFalseColor(p_massremainder[hit_vertices_idx[iv]], plot_v_min, plot_v_max);
                    break;
                case SCMDeformableTerrain::PLOT_ISLAND_ID:
                    mcolor = ChColor(0, 0, 1);
                    if (p_erosion[hit_vertices_idx[iv]] == true)
                        mcolor = ChColor(1, 1, 1);
                    if (p_id_island[hit_vertices_idx[iv]] > 0)
                        mcolor = ChColor::ComputeFalseColor(4.0 + (p_id_island[hit_vertices_idx[iv]] % 8), 0.0, 12.0);
                    if (p_id_island[hit_vertices_idx[iv]] < 0)
                        mcolor = ChColor(0, 0, 0);
                    break;
                case SCMDeformableTerrain::PLOT_IS_TOUCHED:
                    if (p_sigma[hit_vertices_idx[iv]] > 0)
                        mcolor = ChColor(1, 0, 0);
                    else
                        mcolor = ChColor(0, 0, 1);
                    break;
            }
            colors[iv] = {mcolor.R, mcolor.G, mcolor.B};
        }
    } else {
        colors.clear();
    }


    for(int i = 0; i<active_sub_mesh.size(); i++){
        for(int j = 0; j<hit_vertices_idx.size();j++){

            if(hit_vertices_idx[j]<=activeSubMesh_size_buffer[i]){
                if(i==0){
                    m_grid_shape->Update(processed_vertices_hit[j],hit_vertices_idx[j],active_sub_mesh[i]);
                }else{
                    m_grid_shape->Update(processed_vertices_hit[j],hit_vertices_idx[j]-activeSubMesh_size_buffer[i-1]-1,active_sub_mesh[i]);
                }

                if(i==0){
                    m_grid_shape->UpdateColor(colors[j],hit_vertices_idx[j],active_sub_mesh[i]);
                }else{
                    m_grid_shape->UpdateColor(colors[j],hit_vertices_idx[j]-activeSubMesh_size_buffer[i-1]-1,active_sub_mesh[i]);
                }
                
                hit_vertices_idx.erase(hit_vertices_idx.begin()+j);
                original_vertice_hit.erase(original_vertice_hit.begin()+j);
                processed_vertices_hit.erase(processed_vertices_hit.begin()+j);
                colors.erase(colors.begin()+j);
                j = j-1;
            }
        }
    }



    m_timer_visualization.stop();

    // ----------------------------------------------
    // Rectangular Mesh Refinement
    // ----------------------------------------------
    /*
    for(int i = 0; i<hit_vertices_idx.size();i++){
        for(int j = 0; j<active_sub_mesh.size(); j++){
            if(hit_vertices_idx[i]>activeSubMesh_size_buffer[j]){
                m_grid_shape->Refine(vertices[hit_vertices_idx[i]], j+1);
                break;
            }

            if(j==active_sub_mesh.size()-1){
                m_grid_shape->Refine(vertices[hit_vertices_idx[i]], 0);
            }
        }
        
    }
    */
    // ----------------------------------------------
    // ----------------------------------------------


m_timer_total.stop();

}

void SCMDeformableSoilGrid::SetupAuxData(){

    p_vertices_initial = vertices;
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
        //p_level[i] = (vertices[i]).z();
        p_level_initial[i] = p_level[i];
    }
}

void SCMDeformableSoilGrid::GetMesh(){
    m_grid_shape->GetVisMesh(m_trimesh_shape, plane, active_sub_mesh);
}

void SCMDeformableSoilGrid::UpdateFixedPatch(){
    MovingPatchInfo pinfo;
    
    ChVector2<> p_min(+std::numeric_limits<double>::max());
    ChVector2<> p_max(-std::numeric_limits<double>::max());

    // Get current bounding box (AABB) of all collision shapes
    ChVector<> aabb_min;
    ChVector<> aabb_max;
    GetSystem()->GetCollisionSystem()->GetBoundingBox(aabb_min, aabb_max);

    // Loop over all corners of the AABB
    for (int j = 0; j < 8; j++) {
        int ix = j % 2;
        int iy = (j / 2) % 2;
        int iz = (j / 4);

        // AABB corner in absolute frame
        ChVector<> c_abs = aabb_max * ChVector<>(ix, iy, iz) + aabb_min * ChVector<>(1.0 - ix, 1.0 - iy, 1.0 - iz);
        // AABB corner in SCM frame
        ChVector<> c_scm = plane.TransformPointParentToLocal(c_abs);

        // Update AABB of patch projection onto SCM plane
        p_min.x() = std::min(p_min.x(), c_scm.x());
        p_min.y() = std::min(p_min.y(), c_scm.y());
        p_max.x() = std::max(p_max.x(), c_scm.x());
        p_max.y() = std::max(p_max.y(), c_scm.y());

        pinfo.m_max.x()=p_max.x();
        pinfo.m_max.y()=p_max.y();
        pinfo.m_min.x()=p_min.x();
        pinfo.m_min.y()=p_min.y();
        

       m_patches.push_back(pinfo);
    }
}


std::vector<int> FindActiveSubMeshIdx(std::vector<double> x_cut, 
                                    std::vector<double> y_cut, 
                                    std::vector<ChSubGridMeshConnected> subMesh,
                                    std::vector<SCMDeformableSoilGrid::MovingPatchInfo> patches){
    std::vector<int> returnBuffer;


    for(int i = 0; i<subMesh.size();i++){
        
        double sub_cen_x = (subMesh[i].xmin + subMesh[i].xmax) / 2;
        double sub_cen_y = (subMesh[i].ymin + subMesh[i].ymax) / 2;

        for(int j = 0; j<patches.size();j++){
            double patch_cen_x = (patches[j].m_min.x() + patches[j].m_max.x()) / 2;
            double patch_cen_y = (patches[j].m_min.y() + patches[j].m_max.y()) / 2;


            if(patch_cen_x > subMesh[i].xmin && patch_cen_x <subMesh[i].xmax && patch_cen_y > subMesh[i].ymin && patch_cen_y < subMesh[i].ymax)
            {
                returnBuffer.push_back(i);
                break;
            }

            if(sub_cen_x > patches[j].m_min.x() && sub_cen_x <patches[j].m_max.x() && sub_cen_y > patches[j].m_min.y() && sub_cen_y < patches[j].m_max.y())
            {

                returnBuffer.push_back(i);
                break;
            }


            if(patches[j].m_min.x()<subMesh[i].xmax && subMesh[i].xmin<patches[j].m_min.x() && patches[j].m_min.y()<subMesh[i].ymax && patches[j].m_min.y()>subMesh[i].ymin )
            {
                returnBuffer.push_back(i);
                break;
            }


            if(patches[j].m_max.x()<subMesh[i].xmax && subMesh[i].xmin<patches[j].m_max.x() && patches[j].m_min.y()<subMesh[i].ymax && patches[j].m_min.y()>subMesh[i].ymin )
            {

                returnBuffer.push_back(i);
                break;

            }


            if(patches[j].m_max.x()<subMesh[i].xmax && subMesh[i].xmin<patches[j].m_max.x() && patches[j].m_max.y()<subMesh[i].ymax && patches[j].m_max.y()>subMesh[i].ymin )
            {

                returnBuffer.push_back(i);
                break;

            }

            if(patches[j].m_min.x()<subMesh[i].xmax && subMesh[i].xmin<patches[j].m_min.x() && patches[j].m_max.y()<subMesh[i].ymax && patches[j].m_max.y()>subMesh[i].ymin )
            {

                returnBuffer.push_back(i);
                break;

            }

            if(patches[j].m_min.x()<subMesh[i].xmin  && patches[j].m_max.x()>subMesh[i].xmax){
                if(patches[j].m_min.y()<subMesh[i].ymin  && patches[j].m_max.y()>subMesh[i].ymax){

                    returnBuffer.push_back(i);
                    break;

                }
            }

            if(patches[j].m_max.y()>subMesh[i].ymax  && patches[j].m_min.y()<subMesh[i].ymin){
                if(subMesh[i].xmin<patches[j].m_min.x()  && subMesh[i].xmax>patches[j].m_min.x()){

                    returnBuffer.push_back(i);
                    break;

                }
            }

            if(patches[j].m_max.y()>subMesh[i].ymax  && patches[j].m_min.y()<subMesh[i].ymin){
                if(subMesh[i].xmin<patches[j].m_max.x()  && subMesh[i].xmax>patches[j].m_max.x()){

                    returnBuffer.push_back(i);
                    break;
                }
            }


            if(patches[j].m_max.x()>subMesh[i].xmax  && patches[j].m_min.x()<subMesh[i].xmin){
                if(subMesh[i].ymin<patches[j].m_max.y()  && subMesh[i].ymax>patches[j].m_max.y()){

                    returnBuffer.push_back(i);
                    break;

                }
            }

            if(patches[j].m_max.x()>subMesh[i].xmax  && patches[j].m_min.x()<subMesh[i].xmin){
                if(subMesh[i].ymin<patches[j].m_min.y()  && subMesh[i].ymax>patches[j].m_min.y()){

                    returnBuffer.push_back(i);
                    break;

                }
            }
        }
    }

    

    return returnBuffer;
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
