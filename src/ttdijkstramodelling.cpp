/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
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

Dijkstra::Dijkstra( const Graph & graph ) : graph_( graph ) {
    pathMatrix_.resize( graph.size() );
}

void Dijkstra::setGraph( const Graph & graph ) {
    graph_ = graph;
    pathMatrix_.resize( graph.size() );
}

void Dijkstra::setStartNode( int startNode ) {
    distances_.clear();
    root_ = startNode;
    std::priority_queue< distancePair_, std::vector< distancePair_ >, comparePairsClass_< distancePair_ > > priQueue;

    edge_ e( startNode, startNode );
    priQueue.push( distancePair_( 0.0, e ) );
    distancePair_ dummy;

    while ( !priQueue.empty() ) {
        dummy = priQueue.top();

        float distance = dummy.first;
        int node = dummy.second.end;
        priQueue.pop();

        if ( distances_.count( node ) == 0 ) {
            distances_[ node ] = distance;
            pathMatrix_[ node ] = edge_( dummy.second );

            NodeDistMap::iterator start = graph_[ node ].begin();
            NodeDistMap::iterator stop = graph_[ node ].end();

            for (; start != stop; ++start ) {
                e.start = node;
                e.end = (*start).first;
                priQueue.push( distancePair_( distance + (*start).second, e ) );
            }
        }
    }
}

std::vector < int > Dijkstra::shortestPathTo( int node ) const {
    std::vector < int > way;

    int parentNode = -1, endNode = node;

    while ( parentNode != root_ ) {
        parentNode = pathMatrix_[ endNode ].start;
        way.push_back( endNode );
        endNode = pathMatrix_[ endNode ].start;
    }

    way.push_back( root_ );

    std::vector < int > rway( way.size() );
    for ( uint i = 0; i < way.size(); i ++ ) rway[ i ] = way[ way.size() - i - 1 ];

    return rway;
}


//    RVector TravelTimeDijkstraModelling::operator () ( const RVector & slowness, double background ) {
//        return response( slowness, background );
//    }

TravelTimeDijkstraModelling::TravelTimeDijkstraModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose )
: ModellingBase( dataContainer, verbose ), background_( 1e16 ) {

    this->setMesh( mesh );

    this->initJacobian();
}

RVector TravelTimeDijkstraModelling::createDefaultStartModel( ) {
    RVector vec( this->regionManager().parameterCount(), findMedianSlowness() );
    return vec;
}

Graph TravelTimeDijkstraModelling::createGraph( ) {
    Graph meshGraph;

    double dist, oldTime, newTime;
    for ( uint i = 0; i < mesh_->cellCount(); i ++ ) {
        for ( uint j = 0; j < mesh_->cell( i ).nodeCount(); j ++ ) {
                Node *na = &mesh_->cell( i ).node( j );
                Node *nb = &mesh_->cell( i ).node( ( j + 1 )%mesh_->cell( i ).nodeCount() );
                dist = na->pos().distance( nb->pos() );

                oldTime = meshGraph[ na->id() ][ nb->id() ];
                newTime = dist * mesh_->cell( i ).attribute();
                if ( oldTime != 0 ) {
                    newTime = std::min( newTime, oldTime );
                }
                meshGraph[ na->id() ][ nb->id() ] = newTime;
                meshGraph[ nb->id() ][ na->id() ] = newTime;
        }
        if ( mesh_->cell( i ).rtti() == MESH_TETRAHEDRON_RTTI ||
             mesh_->cell( i ).rtti() == MESH_TETRAHEDRON10_RTTI ) {

            Node *na = &mesh_->cell( i ).node( 0 );
            Node *nb = &mesh_->cell( i ).node( 2 );
                dist = na->pos().distance( nb->pos() );
                oldTime = meshGraph[ na->id() ][ nb->id() ];
                newTime = dist * mesh_->cell( i ).attribute();
                meshGraph[ na->id() ][ nb->id() ] = newTime;
                meshGraph[ nb->id() ][ na->id() ] = newTime;

                na = &mesh_->cell( i ).node( 1 );
                nb = &mesh_->cell( i ).node( 3 );
                dist = na->pos().distance( nb->pos() );
                oldTime = meshGraph[ na->id() ][ nb->id() ];
                newTime = dist * mesh_->cell( i ).attribute();
                meshGraph[ na->id() ][ nb->id() ] = newTime;
                meshGraph[ nb->id() ][ na->id() ] = newTime;
        }
    }
    return meshGraph;
}

double TravelTimeDijkstraModelling::findMedianSlowness( ) const {
    return median( getApparentSlowness() );
}

RVector TravelTimeDijkstraModelling::getApparentSlowness( ) const {
    int nData = dataContainer_->size();
    int s = 0, g = 0;
    double edgeLength = 0.0;
    RVector apparentSlowness ( nData );

    for ( int dataIdx = 0; dataIdx < nData; dataIdx ++ ) {
//        s = (*dataContainer_)( dataIdx ).aIdx();
//        g = (*dataContainer_)( dataIdx ).mIdx();
        s = (*dataContainer_)( "s" )[ dataIdx ];
        g = (*dataContainer_)( "g" )[ dataIdx ];
        edgeLength = dataContainer_->sensorPosition( s ).distance( dataContainer_->sensorPosition( g ) );
        apparentSlowness[ dataIdx ] = dataContainer_->get( "t" )[ dataIdx ] / edgeLength;
    }
    return apparentSlowness;
}

RVector TravelTimeDijkstraModelling::calculate( ) {
    dijkstra_.setGraph( createGraph() );

    std::vector < uint > sources( mesh_->findNodesIdxByMarker( -99 ) ); //hack!!!
    int nShots = sources.size();

    std::vector < RVector > dMap( nShots, RVector( nShots ) );

    for ( int shot = 0; shot < nShots; shot ++ ) {
        dijkstra_.setStartNode( sources[ shot ] );

        for ( int i = 0; i < nShots; i ++ ) {
            dMap[ shot ][ i ] = dijkstra_.distance( sources[ i ] );
        }
    }

    int nData = dataContainer_->size();
    int s = 0, g = 0;
    RVector resp( nData );
    for ( int dataIdx = 0; dataIdx < nData; dataIdx ++ ) {
        s = (*dataContainer_)( "s" )[ dataIdx ];
        g = (*dataContainer_)( "g" )[ dataIdx ];

        resp[ dataIdx ] = dMap[ s ][ g ];
    }

    return  resp;
}

RVector TravelTimeDijkstraModelling::response( const RVector & slowness ) {
    if ( background_ < TOLERANCE ) {
        std::cout << "Background: " << background_ << " ->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    this->mapModel( slowness, background_ );
    dijkstra_.setGraph( createGraph() );

    std::vector < uint > sources( mesh_->findNodesIdxByMarker( -99 ) ); //hack!!!
    int nShots = sources.size();

    std::vector < RVector > dMap( nShots, RVector( nShots ) );

    for ( int shot = 0; shot < nShots; shot ++ ) {
        dijkstra_.setStartNode( sources[ shot ] );

        for ( int i = 0; i < nShots; i ++ ) {
            dMap[ shot ][ i ] = dijkstra_.distance( sources[ i ] );
        }
    }

    int nData = dataContainer_->size();
    int s = 0, g = 0;

    RVector resp( nData );

    for ( int dataIdx = 0; dataIdx < nData; dataIdx ++ ) {
        s = (*dataContainer_)( "s" )[dataIdx];
        g = (*dataContainer_)( "g" )[dataIdx];

        resp[ dataIdx ] = dMap[ s ][ g ];
    }

    return  resp;
}

void TravelTimeDijkstraModelling::initJacobian( ){
    if ( jacobian_ && ownJacobian_ ){
        delete jacobian_;
    }
    jacobian_ = new DSparseMapMatrix();
    ownJacobian_ = true;
}


void TravelTimeDijkstraModelling::createJacobian( const RVector & slowness ) {
    DSparseMapMatrix * jacobian = dynamic_cast < DSparseMapMatrix * > ( jacobian_ );
    this->createJacobian( *jacobian, slowness );
}

void TravelTimeDijkstraModelling::createJacobian( DSparseMapMatrix & jacobian, const RVector & slowness ) {
    if ( background_ < TOLERANCE ) {
        std::cout << "Background: " << background_ << " ->" << 1e16 << std::endl;
        background_ = 1e16;
    }

    this->mapModel( slowness, background_ );
        // oder
        // RVector number( slowness.size() );
        // number.fill( x__ );
        // this->mapModel( number, 0.0 );
        // und dann weiter wie vorher nur ohne die -2
    dijkstra_.setGraph( createGraph() );

    std::vector < uint > sources( mesh_->findNodesIdxByMarker( -99 ) ); //hack!!!
    if ( sources.size() < 1 ) {
        //** was passiert hier bei free-elecs oder CEM
        //** (schnellste) Zelle suchen in der Knoten ist und Laufzeit in Zellknoten ausrechnen. Damit weiter dijkstra
        //** oder Hilfsknoten einführen und edges einführen
        THROW_TO_IMPL
    }

    int nShots = sources.size();
    if ( verbose_ ) std::cout << "Found " << nShots << " sources." << std::endl;

    uint nData = dataContainer_->size();
    uint nModel = slowness.size();

    jacobian.clear();
    jacobian.setRows( nData );
    jacobian.setCols( nModel );

    //** for each shot: vector<  way( shot->geoph ) >;
    std::vector < std::vector < std::vector < int > > > wayMatrix( nShots );

    for ( int shot = 0; shot < nShots; shot ++ ) {
        dijkstra_.setStartNode( sources[ shot ] );

        for ( int i = 0; i < nShots; i ++ ) {
            wayMatrix[ shot ].push_back( dijkstra_.shortestPathTo( sources[ i ] ) );
        }
    }

    std::fstream file;
        if ( verbose_ ) openOutFile( "jacobian.way", &file );

        for ( uint dataIdx = 0; dataIdx < nData; dataIdx ++ ) {
            int s = (*dataContainer_)( "s" )[ dataIdx ];
            int g = (*dataContainer_)( "g" )[ dataIdx ];
            std::set < Cell * > neighborCells;

            for ( uint i = 0; i < wayMatrix[ s ][ g ].size()-1; i ++ ) {
                int aId = wayMatrix[ s ][ g ][ i ];
                int bId = wayMatrix[ s ][ g ][ i + 1 ];
                double edgeLength = mesh_->node( aId ).pos().distance( mesh_->node( bId ).pos() );
                double slo = 0.0;

                intersectionSet( neighborCells, mesh_->node( aId ).cellSet(), mesh_->node( bId ).cellSet() );

                if ( verbose_ ) {
                    file << wayMatrix[ s ][ g ][ i ] << " " << wayMatrix[ s ][ g ][ i+1 ] << " " << dijkstra_.distance( wayMatrix[ s ][ g ][ i ] )  << std::endl;
                }

                if ( !neighborCells.empty() ) {
                    double mins = 1e16, dequal = 1e-3;
                    int nfast = 0;
                    /*! first detect cells with minimal slowness */
                    for ( std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++ ) {
                        slo = (*it)->attribute();
                        if ( std::fabs( slo / mins -1 ) < dequal ) nfast++; // check for equality
                        else if ( slo < mins ) {
                            nfast = 1;
                            mins = slo;
                        }
                    }
                    /*! now write edgelength divided by two into jacobian matrix */
                    for ( std::set < Cell * >::iterator it = neighborCells.begin(); it != neighborCells.end(); it ++ ) {
                        int marker = (*it)->marker();
                        if ( nfast > 0 ) {
                            slo = (*it)->attribute();
                            if ( ( slo > 0.0 ) && ( std::fabs( slo / mins - 1 ) < dequal ) ) {
                                if ( marker > (int)nModel - 1 ) {
                                    std::cerr << "Warning! request invalid model cell: " << *(*it) << std::endl;
                                } else {
                                    jacobian[ dataIdx ][ marker ] += edgeLength / nfast; //nur wohin??
                                }
                            }
                        }
                    }
                } else { // neighborCells.empty()
                    std::cerr << WHERE_AM_I << " no neighbor cells found for edge: " << aId << " " << bId << std::endl;
                }
            }
        }
    if ( verbose_ ) file.close( );
}

TTModellingWithOffset::TTModellingWithOffset( Mesh & mesh, DataContainer & dataContainer, bool verbose )
: TravelTimeDijkstraModelling( mesh, dataContainer, verbose ) {

    //! find occuring shots, and map them to indices starting from zero
    shots_ = unique( sort( dataContainer.get("s") ) );
    std::cout << "found " << shots_.size() << " shots." << std::endl;
    for ( size_t i = 0 ; i < shots_.size() ; i++ ) {
        shotMap_.insert( std::pair< int, int >( shots_[ i ], i ) );
    }

    //! create new region containing offsets with special marker

    offsetMesh_ = createMesh1D( shots_.size() );
    for ( size_t i = 0 ; i < offsetMesh_.cellCount() ; i++ ) {
        offsetMesh_.cell( i ).setMarker( NEWREGION );
    }

    regionManager().addRegion( NEWREGION, offsetMesh_ );

    this->initJacobian( );
}

TTModellingWithOffset::~TTModellingWithOffset() { }

RVector TTModellingWithOffset::createDefaultStartModel() {
    return cat( TravelTimeDijkstraModelling::createDefaultStartModel(), RVector( shots_.size() ) );
}

RVector TTModellingWithOffset::response( const RVector & model ) {
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

void TTModellingWithOffset::initJacobian( ){
    if ( jacobian_ && ownJacobian_ ){
        delete jacobian_;
    }
    jacobian_ = new H2SparseMapMatrix();
    ownJacobian_ = true;
}

void TTModellingWithOffset::createJacobian( const RVector & model ){

    H2SparseMapMatrix *jacobian = dynamic_cast < H2SparseMapMatrix* > ( jacobian_ );
    //! extract slowness from model and call old function

    RVector slowness( model, 0, model.size() - shots_.size() );
    RVector offsets( model, model.size() - shots_.size(), model.size() );

    TravelTimeDijkstraModelling::createJacobian( jacobian->H1(), slowness );
    jacobian->H2().setRows( dataContainer_->size() );
    jacobian->H2().setCols( offsets.size() );

    //! set 1 entries for the used shot
    RVector shotpos = dataContainer_->get( "s" ); // shot=C1/A

    for ( size_t i = 0; i < dataContainer_->size(); i++ ) {
        jacobian->H2().setVal( i, shotMap_[ shotpos[ i ] ], 1.0 );
    }
}

} // namespace GIMLI{
