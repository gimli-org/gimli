/***************************************************************************
 *   Copyright (C) 2006-2016 by the resistivity.net development team       *
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

#ifndef _GIMLI_TTDIJKSTRAMODDELING__H
#define _GIMLI_TTDIJKSTRAMODDELING__H

#include "gimli.h"
#include "modellingbase.h"
#include "mesh.h"

namespace GIMLI {

//** sorted vector
typedef std::map< int, double > NodeDistMap;
//** sorted matrix
typedef std::map< int, NodeDistMap > Graph;

/*! Dijkstra's shortest path finding*/
class DLLEXPORT Dijkstra {
public:
    Dijkstra(){}

    Dijkstra(const Graph & graph);

    ~Dijkstra(){}

    void setGraph(const Graph & graph);

    void setStartNode(Index startNode);

    std::vector < Index > shortestPathTo(Index node) const;

    inline double distance(int node) { return distances_[node]; }

    class edge_ : std::pair< int, int > {
    public:
        edge_() : start(0), end(0) {}
        edge_(int a, int b) : start(a), end(b) {}
        int start;
        int end;
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
    int root_;
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

    Graph createGraph();

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
    RVector                 shots_;
    std::map< int, int >    shotMap_;
    Mesh                    offsetMesh_;
};


} //namespace GIMLI

#endif
