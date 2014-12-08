/***************************************************************************
 *   Copyright (C) 2006-2014     by the resistivity.net development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
 *   Thomas Günther thomas@resistivity.net                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "ttdijkstramodelling.h"

#define NEWREGION 33333

#include "blockmatrix.h"
#include "datacontainer.h"
#include "elementmatrix.h"
#include "pos.h"
#include "mesh.h"
#include "meshgenerators.h"
#include "numericbase.h"
#include "regionManager.h"
#include "sparsematrix.h"

#include <vector>
#include <queue>
#include <map>

namespace GIMLI {

Dijkstra::Dijkstra(const Graph & graph) : graph_(graph) {
    pathMatrix_.resize(graph.size());
}

void Dijkstra::setGraph(const Graph & graph) {
    graph_ = graph;
    pathMatrix_.resize(graph.size());
}

void Dijkstra::setStartNode(uint startNode) {
    distances_.clear();
    root_ = startNode;
    std::priority_queue< distancePair_, 
                         std::vector< distancePair_ >,
                         comparePairsClass_< distancePair_ > > priQueue;

    edge_ e(startNode, startNode);
    priQueue.push(distancePair_(0.0, e));
    distancePair_ dummy;
    
    while (!priQueue.empty()) {
        dummy = priQueue.top();

        double distance = dummy.first;
        int node = dummy.second.end;
        priQueue.pop();

        if (distances_.count(node) == 0) {
            distances_[node] = distance;
            if (pathMatrix_.size() <= node){
                throwError(1, WHERE_AM_I + " Warning! Dijkstra graph invalid" );
            }
            pathMatrix_[node] = edge_(dummy.second);

            NodeDistMap::iterator start = graph_[node].begin();
            NodeDistMap::iterator stop = graph_[node].end();

            for (; start != stop; ++start) {
                e.start = node;
                e.end = (*start).first;
                priQueue.push(distancePair_(distance + (*start).second, e));
            }
        }
    }
}

std::vector < uint > Dijkstra::shortestPathTo(uint node) const {
    std::vector < uint > way;

    int parentNode = -1, endNode = node;

    while (parentNode != root_) {
        parentNode = pathMatrix_[endNode].start;
        way.push_back(endNode);
        endNode = pathMatrix_[endNode].start;
    }

    way.push_back(root_);

    std::vector < uint > rway(way.size());
    for (uint i = 0; i < way.size(); i ++) rway[i] = way[way.size() - i - 1];

    return rway;
}


//    RVector TravelTimeDijkstraModelling::operator () (const RVector & slowness, double background) {
//        return response(slowness, background);
//    }

TravelTimeDijkstraModelling::TravelTimeDijkstraModelling(Mesh & mesh,
                                                         DataContainer & dataContainer, 
                                                         bool verbose)
: ModellingBase(dataContainer, verbose), background_(1e16) {

    this->setMesh(mesh);

    this->initJacobian();
}

RVector TravelTimeDijkstraModelling::createDefaultStartModel() {
    RVector vec(this->regionManager().parameterCount(), findMedianSlowness());
    return vec;
}

Graph TravelTimeDijkstraModelling::createGraph() {
    Graph meshGraph;

    double dist, oldTime, newTime;
    
    for (Index i = 0; i < mesh_->cellCount(); i ++) {
        for (Index j = 0; j < mesh_->cell(i).nodeCount(); j ++) {
            Node *na = &mesh_->cell(i).node(j);
            Node *nb = &mesh_->cell(i).node((j + 1)%mesh_->cell(i).nodeCount());
            dist = na->pos().distance(nb->pos());
            
            oldTime = meshGraph[na->id()][nb->id()];
            newTime = dist * mesh_->cell(i).attribute();
            
            if (oldTime != 0) {
                newTime = std::min(newTime, oldTime);
            }
            
            meshGraph[na->id()][nb->id()] = newTime;
            meshGraph[nb->id()][na->id()] = newTime;
        }
        
        if (mesh_->cell(i).rtti() == MESH_TETRAHEDRON_RTTI ||
            mesh_->cell(i).rtti() == MESH_TETRAHEDRON10_RTTI) {

            Node *na = &mesh_->cell(i).node(0);
            Node *nb = &mesh_->cell(i).node(2);
        
            dist = na->pos().distance(nb->pos());
            oldTime = meshGraph[na->id()][nb->id()];
            newTime = dist * mesh_->cell(i).attribute();
        
            meshGraph[na->id()][nb->id()] = newTime;
            meshGraph[nb->id()][na->id()] = newTime;

            na = &mesh_->cell(i).node(1);
            nb = &mesh_->cell(i).node(3);
            
            dist = na->pos().distance(nb->pos());
            oldTime = meshGraph[na->id()][nb->id()];
            newTime = dist * mesh_->cell(i).attribute();
            meshGraph[na->id()][nb->id()] = newTime;
            meshGraph[nb->id()][na->id()] = newTime;
        }
    }
    
    if (meshGraph.size() < mesh_->nodeCount()){
        std::cerr << WHERE_AM_I << 
                " there seems to be unassigned nodes within the mesh. Dijkstra Path will be maybe invalid." 
                 << meshGraph.size() << " < " << mesh_->nodeCount() << std::endl;
        for (Index i = 0; i < mesh_->nodeCount(); i ++){
            if (mesh_->node(i).cellSet().empty()){
                std::cout << mesh_->node(i) << std::endl;
            }
        }
    }
    return meshGraph;
}

double TravelTimeDijkstraModelling::findMedianSlowness() const {
    return median(getApparentSlowness());
}

RVector TravelTimeDijkstraModelling::getApparentSlowness() const {
    uint nData = dataContainer_->size();
    SIndex s = 0, g = 0;
    double edgeLength = 0.0;
    RVector apparentSlowness(nData);

    for (uint dataIdx = 0; dataIdx < nData; dataIdx ++) {
        s = (SIndex)(*dataContainer_)("s")[dataIdx];
        g = (SIndex)(*dataContainer_)("g")[dataIdx];
        if (s == g){
            __MS(WHERE_AM_I + ": shot point equals geophon point. " + 
                  "This lead to an invalid apparent slowness. " + 
                  str(s) + "==" + str(g))
            throwError(1, "Aborting" );
        }
        edgeLength = dataContainer_->sensorPosition(s).distance(dataContainer_->sensorPosition(g));
        apparentSlowness[dataIdx] = dataContainer_->get("t")[dataIdx] / edgeLength;
    }
    return apparentSlowness;
}

RVector TravelTimeDijkstraModelling::createGradientModel(double lBound,
                                                         double uBound){
    if (verbose_) std::cout << "Creating Gradient model ..." << std::endl;
    
    RVector appSlowness(getApparentSlowness());
    
    double smi = median(appSlowness);
    if (smi < lBound) smi = lBound * 1.1;
        
    double sma = max(appSlowness) / 2.0;
    if (uBound > 0.0 && sma > uBound) sma = uBound * 0.9;

    Index nModel = regionManager().parameterCount();

    RVector zmid(nModel);
    Mesh paraDomain(regionManager().paraDomain());
    
    int dim = paraDomain.dim() - 1;

    for (Index i = 0; i < paraDomain.cellCount() ; i++) {
        zmid[i] = paraDomain.cell(i).center()[dim];
    }
    
    double zmi = min(zmid);
    double zma = max(zmid);
        
    RVector gradModel(nModel);

    for (Index i = 0; i < gradModel.size(); i++) {
        gradModel[i] = smi * std::exp((zmid[i] - zmi) / (zma - zmi) * std::log(sma / smi));
    }
    
    return gradModel;
}

void TravelTimeDijkstraModelling::updateMeshDependency_(){
    if (verbose_) std::cout << "... looking for shot and receiver positions." << std::endl;
    RVector shots(unique(sort((*dataContainer_)("s"))));
    
    shotNodeId_.resize(shots.size()) ;
    shotsInv_.clear();
    
    for (uint i = 0; i < shots.size(); i ++){
        shotNodeId_[i] = mesh_->findNearestNode(dataContainer_->sensorPosition((Index)shots[i]));
        if ( mesh_->node(shotNodeId_[i]).cellSet().size() == 0){
        }
        shotsInv_[Index(shots[i])] = i;
    }
     
    RVector receiver(unique(sort((*dataContainer_)("g"))));

    receNodeId_.resize(receiver.size()) ;
    receiInv_.clear();
    
    for (uint i = 0; i < receiver.size(); i ++){
        receNodeId_[i] = mesh_->findNearestNode(dataContainer_->sensorPosition(Index(receiver[i])));
        receiInv_[Index(receiver[i])] = i;
    }
}
   
RVector TravelTimeDijkstraModelling::response(const RVector & slowness) {
    if (background_ < TOLERANCE) {
        std::cout << "Background: " << background_ << "->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    this->mapModel(slowness, background_);
    dijkstra_.setGraph(createGraph());

    uint nShots = shotNodeId_.size();
    uint nRecei = receNodeId_.size();
    RMatrix dMap(nShots, nRecei);

    for (uint shot = 0; shot < nShots; shot ++) {
        dijkstra_.setStartNode(shotNodeId_[shot]);
        for (uint i = 0; i < nRecei; i ++) {
            dMap[shot][i] = dijkstra_.distance(receNodeId_[i]);
        }
    }

    int nData = dataContainer_->size();
    Index s = 0, g = 0;

    RVector resp(nData);

    for (int dataIdx = 0; dataIdx < nData; dataIdx ++) {
        s = shotsInv_[Index((*dataContainer_)("s")[dataIdx])];
        g = receiInv_[Index((*dataContainer_)("g")[dataIdx])];
//         std::cout << s << " " << (*dataContainer_)("s")[dataIdx] << " " 
//                   << g << " " << (*dataContainer_)("g")[dataIdx] << std::endl;
        resp[dataIdx] = dMap[s][g];
    }

    return  resp;
}

void TravelTimeDijkstraModelling::initJacobian(){
    if (jacobian_ && ownJacobian_){
        delete jacobian_;
    }
    jacobian_ = new RSparseMapMatrix();
    ownJacobian_ = true;
}


void TravelTimeDijkstraModelling::createJacobian(const RVector & slowness) {
    RSparseMapMatrix * jacobian = dynamic_cast < RSparseMapMatrix * > (jacobian_);
    this->createJacobian(*jacobian, slowness);
}

void TravelTimeDijkstraModelling::createJacobian(RSparseMapMatrix & jacobian, const RVector & slowness) {
    if (background_ < TOLERANCE) {
        std::cout << "Background: " << background_ << " ->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    this->mapModel(slowness, background_);
    dijkstra_.setGraph(createGraph());

    uint nShots = shotNodeId_.size();
    uint nRecei = receNodeId_.size();
    uint nData = dataContainer_->size();
    uint nModel = slowness.size();

    jacobian.clear();
    jacobian.setRows(nData);
    jacobian.setCols(nModel);

    //** for each shot: vector<  way(shot->geoph) >;
    std::vector < std::vector < std::vector < uint > > > wayMatrix(nShots);

    for (uint shot = 0; shot < nShots; shot ++) {
        dijkstra_.setStartNode(shotNodeId_[shot]);

        for (uint i = 0; i < nRecei; i ++) {
            wayMatrix[shot].push_back(dijkstra_.shortestPathTo(receNodeId_[i]));
        }
    }    
    
    std::fstream file;
        if (verbose_) openOutFile("jacobian.way", &file);

        for (uint dataIdx = 0; dataIdx < nData; dataIdx ++) {
            Index s = shotsInv_[Index((*dataContainer_)("s")[dataIdx])];
            Index g = receiInv_[Index((*dataContainer_)("g")[dataIdx])];
            std::set < Cell * > neighborCells;

            for (uint i = 0; i < wayMatrix[s][g].size()-1; i ++) {
                uint aId = wayMatrix[s][g][i];
                uint bId = wayMatrix[s][g][i + 1];
                double edgeLength = mesh_->node(aId).pos().distance(mesh_->node(bId).pos());
                double slo = 0.0;

                intersectionSet(neighborCells, mesh_->node(aId).cellSet(), mesh_->node(bId).cellSet());

                if (verbose_) {
                    file << wayMatrix[s][g][i] << " " << wayMatrix[s][g][i+1] << " " << dijkstra_.distance(wayMatrix[s][g][i])  << std::endl;
                }

                if (!neighborCells.empty()) {
                    double mins = 1e16, dequal = 1e-3;
                    int nfast = 0;
                    /*! first detect cells with minimal slowness */
                    for (std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++) {
                        slo = (*it)->attribute();
                        if (std::fabs(slo / mins -1) < dequal) nfast++; // check for equality
                        else if (slo < mins) {
                            nfast = 1;
                            mins = slo;
                        }
                    }
                    /*! now write edgelength divided by two into jacobian matrix */
                    for (std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++) {
                        int marker = (*it)->marker();
                        if (nfast > 0) {
                            slo = (*it)->attribute();
                            if ((slo > 0.0) && (std::fabs(slo / mins - 1) < dequal)) {
                                if (marker > (int)nModel - 1) {
                                    std::cerr << "Warning! request invalid model cell: " << *(*it) << std::endl;
                                } else {
                                    jacobian[dataIdx][marker] += edgeLength / nfast; //nur wohin??
                                }
                            }
                        }
                    }
                } else { // neighborCells.empty()
                    std::cerr << WHERE_AM_I << " no neighbor cells found for edge: " << aId << " " << bId << std::endl;
                }
            }
        }
    if (verbose_) file.close();
}

TTModellingWithOffset::TTModellingWithOffset(Mesh & mesh, DataContainer & dataContainer, bool verbose)
: TravelTimeDijkstraModelling(mesh, dataContainer, verbose) {

    //! find occuring shots, and map them to indices starting from zero
    shots_ = unique(sort(dataContainer.get("s")));
    std::cout << "found " << shots_.size() << " shots." << std::endl;
    for (Index i = 0 ; i < shots_.size() ; i++) {
        shotMap_.insert(std::pair< int, int >((Index)shots_[i], i));
    }

    //! create new region containing offsets with special marker

    offsetMesh_ = createMesh1D(shots_.size());
    for (size_t i = 0 ; i < offsetMesh_.cellCount() ; i++) {
        offsetMesh_.cell(i).setMarker(NEWREGION);
    }

    regionManager().addRegion(NEWREGION, offsetMesh_);

    this->initJacobian();
}

TTModellingWithOffset::~TTModellingWithOffset() { }

RVector TTModellingWithOffset::createDefaultStartModel() {
    return cat(TravelTimeDijkstraModelling::createDefaultStartModel(), RVector(shots_.size()));
}

RVector TTModellingWithOffset::response(const RVector & model) {
    //! extract slowness from model and call old function

    RVector slowness(model, 0, model.size() - shots_.size());
    RVector offsets(model, model.size() - shots_.size(), model.size());
    RVector resp = TravelTimeDijkstraModelling::response(slowness); //! normal response
    RVector shotpos = dataContainer_->get("s");

    for (size_t i = 0; i < resp.size() ; i++){
        resp[i] += offsets[shotMap_[Index(shotpos[i])]];
    }

    return resp;
}

void TTModellingWithOffset::initJacobian(){
    if (jacobian_ && ownJacobian_){
        delete jacobian_;
    }
    jacobian_ = new H2SparseMapMatrix();
    ownJacobian_ = true;
}

void TTModellingWithOffset::createJacobian(const RVector & model){

    H2SparseMapMatrix *jacobian = dynamic_cast < H2SparseMapMatrix* > (jacobian_);
    //! extract slowness from model and call old function

    RVector slowness(model, 0, model.size() - shots_.size());
    RVector offsets(model, model.size() - shots_.size(), model.size());

    TravelTimeDijkstraModelling::createJacobian(jacobian->H1(), slowness);
    jacobian->H2().setRows(dataContainer_->size());
    jacobian->H2().setCols(offsets.size());

    //! set 1 entries for the used shot
    RVector shotpos = dataContainer_->get("s"); // shot=C1/A

    for (size_t i = 0; i < dataContainer_->size(); i++) {
        jacobian->H2().setVal(i, shotMap_[Index(shotpos[i])], 1.0);
    }
}

} // namespace GIMLI{
