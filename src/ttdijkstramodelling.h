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

#ifndef _GIMLI_TTDIJKSTRAMODDELING__H
#define _GIMLI_TTDIJKSTRAMODDELING__H

#include "gimli.h"
#include "modellingbase.h"
#include "mesh.h"

namespace GIMLI {

typedef std::map< Index, double > NodeDistMap;

typedef std::map< Index, NodeDistMap > Graph;

/*! Dijkstra's shortest path finding*/
class DLLEXPORT Dijkstra {
public:
    Dijkstra(){}

    Dijkstra(const Graph & graph);

    ~Dijkstra(){}

    void setGraph(const Graph & graph);

    void setStartNode(Index startNode);

    IndexArray shortestPathTo(Index node) const;

    inline double distance(Index node) { return distances_[node]; }

    RVector distances() const;

    class edge_ : std::pair< Index, Index > {
    public:
        edge_() : start(0), end(0) {}
        edge_(Index a, Index b) : start(a), end(b) {}
        Index start;
        Index end;
    };

    /*! Definition for the priority queue */
    class distancePair_ : std::pair< float, edge_ > { // weigth, vertex;
    public:
        distancePair_() : first(0.0) {}
        distancePair_(double f, edge_ & s) : first(f), second(s) {}

        double first;
        edge_ second;
    };

    template < class T > class comparePairsClass_ : public std::binary_function< T, T, T > {
    public:
        bool operator() (const T & lhs, const T & rhs) {
            return lhs.first > rhs.first;
        }
    };

protected:
    std::vector < edge_ > pathMatrix_;
    NodeDistMap distances_;
    Graph graph_;
    Index root_;
};

//! Modelling class for travel time problems using the Dijkstra algorithm
/*! TravelTimeDijkstraModelling(mesh, datacontainer) */
class DLLEXPORT TravelTimeDijkstraModelling : public ModellingBase {
public:
    TravelTimeDijkstraModelling(bool verbose=false);

    TravelTimeDijkstraModelling(Mesh & mesh,
                                DataContainer & dataContainer,
                                bool verbose=false);

    virtual ~TravelTimeDijkstraModelling() { }

    virtual RVector createDefaultStartModel();

    RVector createGradientModel(double lBound, double uBound);

    /*! Interface. Calculate response */
    virtual RVector response(const RVector & slowness);

    /*! Interface. */
    virtual void createJacobian(const RVector & slowness);

    /*! Interface. */
    virtual void initJacobian();

    Graph createGraph(const RVector & slownessPerCell) const;

//     RVector calculate();

    double findMedianSlowness() const;

    RVector getApparentSlowness() const;

    void createJacobian(RSparseMapMatrix & jacobian, const RVector & slowness);

protected:

    /*! Automatically looking for shot and receiver points if the mesh is changed. */
    virtual void updateMeshDependency_();

    Dijkstra dijkstra_;
    double background_;

    /*! Nearest nodes for the current mesh for all shot points.*/
    std::vector < Index > shotNodeId_;

    /*! Map shot id to sequential shot node number of shotNodeId_ */
    std::map< Index, Index > shotsInv_;

    /*! Nearest nodes for the current mesh for all receiver points.*/
    std::vector < Index > receNodeId_;

    /*! Map receiver id to sequential receiver node number of receNodeId_ */
    std::map< Index, Index > receiInv_;

};

/*! New Class derived from standard travel time modelling */
class DLLEXPORT TTModellingWithOffset: public TravelTimeDijkstraModelling{
public:
    TTModellingWithOffset(Mesh & mesh, DataContainer & dataContainer, bool verbose);

    virtual ~TTModellingWithOffset();

    virtual RVector createDefaultStartModel();

    virtual RVector response(const RVector & model);

    void initJacobian();

    virtual void createJacobian(const RVector & slowness);

    size_t nShots(){ return shots_.size(); }

protected:
    RVector                  shots_;
    std::map< Index, Index > shotMap_;
    Mesh                     offsetMesh_;
};


} //namespace GIMLI

#endif
