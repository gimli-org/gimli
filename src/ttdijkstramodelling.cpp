/******************************************************************************
 *   Copyright (C) 2006-2018 by the GIMLi development team                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *   Thomas Günther thomas@resistivity.net                                    *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

#include "ttdijkstramodelling.h"

#define NEWREGION 33333

#include "blockmatrix.h"
#include "datacontainer.h"
#include "elementmatrix.h"
#include "pos.h"
#include "mesh.h"
#include "meshgenerators.h"
#include "meshentities.h"
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

double Dijkstra::distance(Index node) { 
    return distances_[node].time(); 
}

RVector Dijkstra::distances(bool withSecNodes) const {
    RVector ret(0);
    for (auto const & it: distances_){
        ret.push_back(it.second.time());
    }
    return ret;
}

void Dijkstra::setGraph(const Graph & graph) {
    graph_ = graph;
    pathMatrix_.clear();
    pathMatrix_.resize(graph.size());
}

void Dijkstra::setStartNode(Index startNode) {
    distances_.clear();
    root_ = startNode;
    std::priority_queue< DistancePair_,
                         std::vector< DistancePair_ >,
                         ComparePairsClass_< DistancePair_ > > priQueue;

    Edge_ e(startNode, startNode);
    priQueue.push(DistancePair_(0.0, e));
    DistancePair_ dummy;

    while (!priQueue.empty()) {
        dummy = priQueue.top();

        double distance = dummy.first;
        Index node = dummy.second.end;
        priQueue.pop();

        if (distances_.count(node) == 0) {
            //distances_[node] = distance;
            
            distances_[node] = GraphDistInfo(distance, 0);

            if ((Index)pathMatrix_.size() <= node){
                std::cout << "startNodeID:" << startNode << " NodeID:" << node << std::endl;
                throwError(1, WHERE_AM_I + " Warning! Dijkstra graph invalid" );
            }
            pathMatrix_[node] = Edge_(dummy.second);

            NodeDistMap::iterator start = graph_[node].begin();
            NodeDistMap::iterator stop = graph_[node].end();

            for (; start != stop; ++start) {
                e.start = node;
                e.end = (*start).first;
                priQueue.push(DistancePair_(distance + (*start).second.time(), e));
            }
        }
    }
}

IndexArray Dijkstra::shortestPathTo(Index node) const {
    IndexArray way;

    Index parentNode = -1, endNode = node;

    while (parentNode != root_) {
        parentNode = pathMatrix_[endNode].start;
        way.push_back(endNode);
        endNode = pathMatrix_[endNode].start;
    }

    way.push_back(root_);

    IndexArray rway(way.size());
    for (Index i = 0; i < way.size(); i ++) rway[i] = way[way.size() - i - 1];

    return rway;
}


//    RVector TravelTimeDijkstraModelling::operator () (const RVector & slowness, double background) {
//        return response(slowness, background);
//    }

TravelTimeDijkstraModelling::TravelTimeDijkstraModelling(bool verbose)
    : ModellingBase(verbose), background_(1e16){
    this->initJacobian();
}

TravelTimeDijkstraModelling::TravelTimeDijkstraModelling(Mesh & mesh,
                                                         DataContainer & dataContainer,
                                                         bool verbose)
    : ModellingBase(dataContainer, verbose), background_(1e16) {

    this->setMesh(mesh);
    this->initJacobian();
}

RVector TravelTimeDijkstraModelling::createDefaultStartModel() {
    return RVector(this->regionManager().parameterCount(), findMedianSlowness());
}

bool V_ = false;

void fillGraph_(Graph & graph, const Node & a, const Node & b, double slowness, SIndex leftID){
    if (a.id() == b.id()) return;

    double dist = a.pos().distance(b.pos());

    // ensure connection between 3d boundaries
    dist = max(1e-8, dist);

    double newTime = dist * slowness;
    double oldTime = graph[a.id()][b.id()].time();
    
    // if (V_){
    //     __MS("a:" << a.id() << " b:"  << b.id() << " L:" << leftID << " t:" << " " << newTime << " " << oldTime)
    // }

    if (oldTime > 0.0) {
        newTime = std::min(newTime, oldTime);

        // way pair already exist so set time to min and add leftID

        NodeDistMap::iterator ita(graph[a.id()].find(b.id()));
        ita->second.cellIDs().insert(leftID);
        ita->second.setTime(newTime);
        
        NodeDistMap::iterator itb(graph[b.id()].find(a.id()));
        itb->second.cellIDs().insert(leftID);
        itb->second.setTime(newTime);
        
    } else {
        // first time fill
        graph[a.id()][b.id()] = GraphDistInfo(newTime, dist, leftID);   
        graph[b.id()][a.id()] = GraphDistInfo(newTime, dist, leftID);    
    }
}

void fillGraph_(Graph & graph, Cell & c, double slowness){

    std::vector< Node * > ni(c.nodes());

    for (Index i(0); i < c.boundaryCount(); i++){
        Boundary *b = c.boundary(i);
        if (b){
            for (auto & n : b->secondaryNodes()){
                ni.push_back(n);
            }
        } else {
            log(Critical, "No boundary found.");
        }
    }

    for (auto & n : c.secondaryNodes()){
        ni.push_back(n);
    }

    for (Index j = 0; j < ni.size()-1; j ++) {
        for (Index k = j + 1; k < ni.size(); k ++) {
            fillGraph_(graph, *ni[j], *ni[k], slowness, c.id());
        }
    }
}

Graph TravelTimeDijkstraModelling::createGraph(const RVector & slownessPerCell) const {
    Graph graph;
    mesh_->createNeighbourInfos();

    for (Index i = 0; i < mesh_->cellCount(); i ++) {
        Cell & c = mesh_->cell(i);
        fillGraph_(graph, c, slownessPerCell[c.id()]);
    }

    if (graph.size() < mesh_->nodeCount()){
        std::cerr << WHERE_AM_I <<
                " there seems to be unassigned nodes within the mesh. Dijkstra Path will be maybe invalid."
                 << graph.size() << " < " << mesh_->nodeCount() << std::endl;
    }
    return graph;
}

double TravelTimeDijkstraModelling::findMedianSlowness() const {
    return median(getApparentSlowness());
}

RVector TravelTimeDijkstraModelling::getApparentSlowness() const {
    if (!dataContainer_) return 0.0;

    Index nData = dataContainer_->size();
    SIndex s = 0, g = 0;
    double edgeLength = 0.0;
    RVector apparentSlowness(nData);

    for (Index dataIdx = 0; dataIdx < nData; dataIdx ++) {
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

    if (!dataContainer_){
        throwError(1, "We have no dataContainer defined");
    }
    RVector shots(unique(sort((*dataContainer_)("s"))));

    if (shots.size() == 0){
        throwError(1, "There are no shot positions in the dataContainer.");
    }
    shotNodeId_.resize(shots.size()) ;
    shotsInv_.clear();

    if (shots[0] < 0){
        throwError(1, "There are shots index lower then 0.");
    }

    for (Index i = 0; i < shots.size(); i ++){
        shotNodeId_[i] = mesh_->findNearestNode(dataContainer_->sensorPosition((Index)shots[i]));
        if (mesh_->node(shotNodeId_[i]).cellSet().size() == 0){
            __MS("no cells found")
        }
        shotsInv_[Index(shots[i])] = i;
    }

    RVector receiver(unique(sort((*dataContainer_)("g"))));

    receNodeId_.resize(receiver.size()) ;
    receiInv_.clear();

    for (Index i = 0; i < receiver.size(); i ++){
        receNodeId_[i] = mesh_->findNearestNode(dataContainer_->sensorPosition(Index(receiver[i])));
        receiInv_[Index(receiver[i])] = i;
    }
}

RVector TravelTimeDijkstraModelling::response(const RVector & slowness) {
    if (background_ < TOLERANCE) {
        std::cout << "Background: " << background_ << "->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    RVector slowPerCell(this->createMappedModel(slowness, background_));
    dijkstra_.setGraph(createGraph(slowPerCell));

    // this->mapModel(slowness, background_);
    // dijkstra_.setGraph(createGraph(mesh_->cellAttributes()));

    Index nShots = shotNodeId_.size();
    Index nRecei = receNodeId_.size();
    RMatrix dMap(nShots, nRecei);

    for (Index shot = 0; shot < nShots; shot ++) {
        dijkstra_.setStartNode(shotNodeId_[shot]);
        for (Index i = 0; i < nRecei; i ++) {
            dMap[shot][i] = dijkstra_.distance(receNodeId_[i]);
        }
    }

    Index nData = dataContainer_->size();
    Index s = 0, g = 0;

    RVector resp(nData);

    for (Index dataIdx = 0; dataIdx < nData; dataIdx ++) {
        s = shotsInv_[Index((*dataContainer_)("s")[dataIdx])];
        g = receiInv_[Index((*dataContainer_)("g")[dataIdx])];
//         if (dataIdx < 10 ) std::cout << s << " " << (*dataContainer_)("s")[dataIdx] << " "
//                    << g << " " << (*dataContainer_)("g")[dataIdx] << " " << dMap[s][g] << std::endl;
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

const IndexArray & TravelTimeDijkstraModelling::way(Index sht, Index rec) const{
    ASSERT_SIZE(wayMatrix_, sht)
    ASSERT_SIZE(wayMatrix_[sht], rec)
    return wayMatrix_[sht][rec];
}

void TravelTimeDijkstraModelling::createJacobian(const RVector & slowness) {
    RSparseMapMatrix * jacobian = dynamic_cast < RSparseMapMatrix * > (jacobian_);
    this->createJacobian(*jacobian, slowness);
}

void TravelTimeDijkstraModelling::createJacobian(RSparseMapMatrix & jacobian,
                                                 const RVector & slowness) {
    if (background_ < TOLERANCE) {
        std::cout << "Background: " << background_ << " ->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    RVector slowPerCell(this->createMappedModel(slowness, background_));
    dijkstra_.setGraph(createGraph(slowPerCell));

    Index nShots = shotNodeId_.size();
    Index nRecei = receNodeId_.size();
    Index nData = dataContainer_->size();
    Index nModel = slowness.size();

    jacobian.clear();
    jacobian.setRows(nData);
    jacobian.setCols(nModel);

    //** for each shot: vector<  way(shot->geoph) >;
    wayMatrix_.clear();
    wayMatrix_.resize(nShots);

    for (Index shot = 0; shot < nShots; shot ++) {
        dijkstra_.setStartNode(shotNodeId_[shot]);

        for (Index i = 0; i < nRecei; i ++) {
            wayMatrix_[shot].push_back(dijkstra_.shortestPathTo(receNodeId_[i]));
        }
    }

    for (Index dataIdx = 0; dataIdx < nData; dataIdx ++) {
        Index s = shotsInv_[Index((*dataContainer_)("s")[dataIdx])];
        Index g = receiInv_[Index((*dataContainer_)("g")[dataIdx])];

        std::set < Cell * > neighborCells;

        for (Index i = 0; i < wayMatrix_[s][g].size()-1; i ++) {
            neighborCells.clear();

            Index aId = wayMatrix_[s][g][i];
            Index bId = wayMatrix_[s][g][i + 1];
       
            const GraphDistInfo & way = dijkstra_.graphInfo(aId, bId);
            
            double edgeLength = way.dist();
            //double edgeLength = mesh_->node(aId).pos().distance(mesh_->node(bId).pos());
            double slo = 0.0;

            double minSlow = 9e99;

            for (const auto &iCD : way.cellIDs()){
                minSlow = min(minSlow, slowPerCell[iCD]);
            }
            
            for (const auto &iCD : way.cellIDs()){
                if (std::fabs(slowPerCell[iCD] - minSlow) < 1e-4){
                    Cell *c = & mesh_->cell(iCD);
                    neighborCells.insert(c);
                }
            }

            for (const auto &c : neighborCells){
                jacobian[dataIdx][c->marker()] += edgeLength / neighborCells.size();
            } 

//             continue;


//             intersectionSet(neighborCells, mesh_->node(aId).cellSet(), mesh_->node(bId).cellSet());
//             // __MS(aId << " " << bId)

//             if (!neighborCells.empty()) {
//                 double mins = 1e16, dEqual = 1e-3;
//                 int nfast = 0;
//                 /*! first detect cells with minimal slowness */
//                 for (std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++) {
//                     slo = slowPerCell[(*it)->id()];
//                     if (std::fabs(slo / mins -1) < dEqual) nfast++; // check for equality
//                     else if (slo < mins) {
//                         nfast = 1;
//                         mins = slo;
//                     }
//                 }
                
//                 // __MS(aId << " " << bId << "" << slo)

//                 /*! now write edgelength divided by two into jacobian matrix */
//                 for (std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++) {
//                     int marker = (*it)->marker();
//                     if (nfast > 0) {
//                         slo = slowPerCell[(*it)->id()];
//                         if ((slo > 0.0) && (std::fabs(slo / mins - 1) < dEqual)) {
//                             if (marker > (int)nModel - 1) {
//                                 std::cerr << "Warning! request invalid model cell: " << *(*it) << std::endl;
//                             } else {

//                                 if (marker <= MARKER_FIXEDVALUE_REGION){
//                                     // neighbor is fixed region
// //                                         SIndex regionMarker = -(marker - MARKER_FIXEDVALUE_REGION);
// //                                         double val = regionManager_->region(regionMarker)->fixValue();
//                                 } else {
//                                     __MS(aId << " " << bId << " " << marker << " "<< edgeLength / nfast)
//                                     jacobian[dataIdx][marker] += edgeLength / nfast; //nur wohin?? CA nur wohin was??
//                                 }
//                             }
//                         }
//                     }
//                 }
//             } else { // neighborCells.empty()
//                 std::cerr << WHERE_AM_I << " no neighbor cells found for edge: " << aId << " " << bId << std::endl;
//             }
        }
    }
}

TTModellingWithOffset::TTModellingWithOffset(Mesh & mesh, DataContainer & dataContainer, bool verbose)
: TravelTimeDijkstraModelling(mesh, dataContainer, verbose) {

    //! find occuring shots, and map them to indices starting from zero
    shots_ = unique(sort(dataContainer.get("s")));
    std::cout << "found " << shots_.size() << " shots." << std::endl;
    for (Index i = 0 ; i < shots_.size() ; i++) {
        shotMap_.insert(std::pair< Index, Index >((Index)shots_[i], i));
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
