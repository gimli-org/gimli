/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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
#include "regionmanager.h"
#include "dataContainer.h"
#include "blockmatrix.h"
#include "meshgenerators.h"
#define NEWREGION 33333

namespace GIMLI {
    class Mesh;
    class DataMap;

//** sorted vector
typedef std::map< int, double > NodeDistMap;
//** sorted matrix
typedef std::map< int, NodeDistMap > Graph;

/*! Dijkstra's shortest path finding*/
class DLLEXPORT Dijkstra {
public:
    Dijkstra( ){}

    Dijkstra( const Graph & graph );

    ~Dijkstra( ){}

    void setGraph( const Graph & graph );

    void setStartNode( int startNode );

    std::vector < int > shortestPathTo( int node ) const;

    inline double distance( int node ) { return distances_[ node ]; }

    class edge_ : std::pair< int, int > {
    public:
        edge_( ) : start( 0 ), end( 0 ) {}
        edge_( int a, int b ) : start( a ), end( b ) {}
        int start;
        int end;
    };

    /*! Definition for the priority queue */
    class distancePair_ : std::pair< float, edge_ > { // weigth, vertex;
    public:
        distancePair_() : first( 0.0 ) {}
        distancePair_( double f, edge_ & s ) : first( f ), second( s ) {}

        double first;
        edge_ second;
    };

    template < class T > class comparePairsClass_ : public std::binary_function< T, T, T > {
    public:
        bool operator() ( const T & lhs, const T & rhs) {
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
/*! TravelTimeDijkstraModelling( mesh, datacontainer ) */
class DLLEXPORT TravelTimeDijkstraModelling : public ModellingBase {
public:
    TravelTimeDijkstraModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose = false )
        : ModellingBase( dataContainer, verbose ), background_( 1e16 ) {
        setMesh( mesh );
    }

    virtual ~TravelTimeDijkstraModelling() { }

    RVector createDefaultStartModel( );

    /*! Calculate response */
    virtual RVector response( const RVector & slowness );

//    RVector response( const RVector & slowness, double background = 1e16 );
//    virtual RVector operator () ( const RVector & slowness ); //! shortcut for f( model )

    Graph createGraph( );

    RVector calculate( );

    double findMedianSlowness( ) const;

    RVector getApparentSlowness( ) const;

    void createJacobian( DSparseMapMatrix & jacobian, const RVector & slowness );

protected:
    Dijkstra dijkstra_;
    double background_;
};

/*! New Class derived from standard travel time modelling */
class TTModellingWithOffset: public TravelTimeDijkstraModelling{
public:
    TTModellingWithOffset( Mesh & mesh, DataContainer & dataContainer, bool verbose) :
        TravelTimeDijkstraModelling( mesh, dataContainer, verbose ) {
        //! find occuring shots, and map them to indices starting from zero
        shots_ = unique( sort( dataContainer.get("s") ) );
        std::cout << "found " << shots_.size() << " shots." << std::endl;
        for ( size_t i = 0 ; i < shots_.size() ; i++ ) {
            shotMap_.insert( std::pair< int, int >( shots_[ i ], i ) );
        }
        //! create new region containing offsets with special marker
        offsetMesh_ = createMesh1D( shots_.size() );
        for ( size_t i = 0 ; i < offsetMesh_.cellCount() ; i++ ) offsetMesh_.cell( i ).setMarker( NEWREGION );
        regionManager().addRegion( NEWREGION, offsetMesh_ );
    }

    virtual ~TTModellingWithOffset() { }

    RVector createDefaultStartModel() { 
        return cat( TravelTimeDijkstraModelling::createDefaultStartModel(), RVector( shots_.size() ) );
    }
    
    RVector response( const RVector & model ) {
        //! extract slowness from model and call old function
        RVector slowness( model, 0, model.size() - shots_.size() );
        RVector offsets( model, model.size() - shots_.size(), model.size() ); 
        RVector resp = TravelTimeDijkstraModelling::response( slowness ); //! normal response
        RVector shotpos = dataContainer_->get( "s" );
        for ( size_t i = 0; i < resp.size() ; i++ ){
            resp[ i ] += offsets[ shotMap_[ shotpos[ i ] ] ];
        }
        return resp;
    }

    void createJacobian( H2SparseMapMatrix1 & jacobian, const RVector & model ) { 
        //! extract slowness from model and call old function
        RVector slowness( model, 0, model.size() - shots_.size() );
        RVector offsets( model, model.size() - shots_.size(), model.size() );         
        TravelTimeDijkstraModelling::createJacobian( *jacobian.H1(), slowness );
        jacobian.H2()->setRows( dataContainer_->size() );
        jacobian.H2()->setCols( offsets.size() );
        //! set 1 entries for the used shot
        RVector shotpos = dataContainer_->get( "s" ); // shot=C1/A
        for ( size_t i = 0; i < dataContainer_->size(); i++ ) {
            jacobian.H2()->setVal( i, shotMap_[ shotpos[ i ] ], 1.0 ); 
        }
    }
    
    size_t nShots(){ return shots_.size(); }
protected:
    RVector shots_;
    std::map< int, int > shotMap_;
    Mesh offsetMesh_;
};


} //namespace GIMLI

#endif
