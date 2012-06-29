/***************************************************************************
 *   Copyright (C) 2008-2011 by the resistivity.net development team       *
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

#include "meshgenerators.h"
#include "mesh.h"
#include "meshentities.h"
#include "triangleWrapper.h"

namespace GIMLI{

Mesh createMesh1D( uint nCells, uint nClones ){
    RVector pos( nCells * nClones + 1);

    std::generate( pos.begin(), pos.end(), IncrementSequence< double >( 0.0 ) );
    Mesh mesh( createMesh1D( pos ) );
    for ( uint i = 0; i < nClones; i ++ ) {
        for ( uint j = 0; j < nCells; j ++ ) {
            mesh.cell( ( i * nCells ) + j ).setMarker( i );
        }
    }

  return mesh;
}

Mesh createMesh1D( const std::vector < double > & x ){
    return createMesh1D( RVector( x ) );
}

Mesh createMesh1D( const RVector & x ){
    Mesh mesh( 1 );
    mesh.create1DGrid( x );
    return mesh;
}

Mesh createMesh1DBlock( uint nLayers, uint nProperties ){
    uint nPar = nLayers * ( nProperties + 1 );
    RVector pos( nPar );
    std::generate( pos.begin(), pos.end(), IncrementSequence< double>( 0.0 ) );
    Mesh mesh( createMesh1D( pos ) );
    /*! Thicknesses have marker 0 */
    for ( uint i = 0; i < nLayers - 1; i++ ) mesh.cell( i ).setMarker( 0 );
    /*! Properties have markers 1,2,... */
    for ( uint i = 0; i < nProperties; i ++ ) {
        for ( uint j = 0; j < nLayers; j ++ ) {
            mesh.cell( ( i + 1 ) * nLayers + j -1 ).setMarker( i + 1 );
        }
    }
    return mesh;
}

Mesh createMesh2D( uint xDim, uint yDim, int markerType ){
    RVector x( xDim + 1 ); std::generate( x.begin(), x.end(), IncrementSequence< double >( 0.0 ) );
    RVector y( yDim + 1 ); std::generate( y.begin(), y.end(), IncrementSequence< double >( 0.0 ) );
    return createMesh2D( x, y, markerType );
}

Mesh createMesh2D( const RVector & x, const RVector & y, int markerType ){
    Mesh mesh( 2 );
    mesh.create2DGrid( x, y, markerType );
    for (uint i = 0; i < mesh.boundaryCount(); i ++ ){
        if ( !mesh.boundary( i ).leftCell() || !mesh.boundary( i ).rightCell() ){
            mesh.boundary( i ).setMarker( 1 );
        }
    }
    return mesh;
}

Mesh createMesh3D( uint xDim, uint yDim, uint zDim, int markerType ){
    RVector x( xDim + 1 ); std::generate( x.begin(), x.end(), IncrementSequence< double >( 0.0 ) );
    RVector y( yDim + 1 ); std::generate( y.begin(), y.end(), IncrementSequence< double >( 0.0 ) );
    RVector z( zDim + 1 ); std::generate( z.begin(), z.end(), IncrementSequence< double >( 0.0 ) );
    return createMesh3D( x, y, z, markerType );
}

Mesh createMesh3D( const RVector & x, const RVector & y, const RVector & z, int markerType ){
    Mesh mesh( 3 );
    mesh.create3DGrid( x, y, z, markerType );
    int marker = 1;
    for (uint i = 0; i < mesh.boundaryCount(); i ++ ){
        if ( !mesh.boundary( i ).leftCell() || !mesh.boundary( i ).rightCell() ){
            mesh.boundary( i ).setMarker( marker );
        }
    }
    return mesh;
}

Mesh createMesh3D( const Mesh mesh, const RVector & z, int markerType ){
    Mesh mesh3( 3 );
    
    if ( z.size() < 2 ){
        std::cout << "Warning!: " << WHERE_AM_I << "extrusion vector size need z be greater than 1" << std::endl;
    }
        
    for ( Index iz = 0; iz < z.size(); iz ++ ){
        for ( Index ic = 0; ic < mesh.nodeCount(); ic ++ ){
            mesh3.createNode( mesh.node( ic ).pos() + RVector3( 0.0, 0.0, z[ iz ] ) );
        }
    }
    
    std::vector < Node * > nodes( 6 );
    for ( Index iz = 1; iz < z.size(); iz ++ ){
        for ( Index ic = 0; ic < mesh.cellCount(); ic ++ ){
            // "check for triangle here"
            
            nodes[ 0 ] = & mesh3.node( (iz-1) * mesh.nodeCount() + mesh.cell( ic ).node( 0 ).id() );
            nodes[ 1 ] = & mesh3.node( (iz-1) * mesh.nodeCount() + mesh.cell( ic ).node( 1 ).id() );
            nodes[ 2 ] = & mesh3.node( (iz-1) * mesh.nodeCount() + mesh.cell( ic ).node( 2 ).id() );
            nodes[ 3 ] = & mesh3.node( (iz) * mesh.nodeCount() + mesh.cell( ic ).node( 0 ).id() );
            nodes[ 4 ] = & mesh3.node( (iz) * mesh.nodeCount() + mesh.cell( ic ).node( 1 ).id() );
            nodes[ 5 ] = & mesh3.node( (iz) * mesh.nodeCount() + mesh.cell( ic ).node( 2 ).id() );
            mesh3.createCell( nodes, markerType );
        }
    }
    
    
    return mesh3;
}

bool addTriangleBoundary( Mesh & mesh, double xBoundary, double yBoundary, int cellMarker, 
                          bool save ){

    int boundMarker = -5;
    
    mesh.createNeighbourInfos( true );
    Boundary *b = NULL;
    for ( std::vector < Boundary * >::const_iterator it = mesh.boundaries().begin();
        it != mesh.boundaries().end(); it ++ ){
        b = (*it);
        if ( ! b->leftCell() || ! b->rightCell() ) {
            if ( b->marker() != MARKER_BOUND_HOMOGEN_NEUMANN ){
                b->setMarker( boundMarker );
            }
        }
    }
    
    std::vector < Boundary * > markedBoundaries( mesh.findBoundaryByMarker( boundMarker ) );
    Boundary *start( markedBoundaries.front() );
    
    std::list < Node * > boundNodes;
    
    boundNodes.push_back( & start->node( 0 ) );
    Boundary * thisBoundary = start;
    
    while ( thisBoundary != NULL ){
        Node * thisNode = boundNodes.back();
        Boundary * nextBoundary = NULL;
    
        for ( std::set < Boundary * >::iterator it = thisNode->boundSet().begin();
                it != thisNode->boundSet().end(); it++ ){
        
            if ( (*it)->marker() == boundMarker && (*it) != thisBoundary ){
                nextBoundary = (*it);
                break;
            }
        }
        if ( nextBoundary ){
            Node * nextNode = NULL; 
            if ( &nextBoundary->node( 0 ) != thisNode ){
                nextNode = &nextBoundary->node( 0 );
            } else {
                nextNode = &nextBoundary->node( 1 );
            }
            //** Check closed polygones here 
            if ( find( boundNodes.begin(), boundNodes.end(), nextNode ) == boundNodes.end() ) {
                boundNodes.push_back( nextNode );
            } else {
                break;
            }
        }
        thisBoundary = nextBoundary;
    }
    boundNodes.push_front( & start->node( 1 ) );
    thisBoundary = start;
    
    while ( thisBoundary != NULL ){
        Node * thisNode = boundNodes.front();
        Boundary * nextBoundary = NULL;
    
        for ( std::set < Boundary * >::iterator it = thisNode->boundSet().begin();
                it != thisNode->boundSet().end(); it++ ){
        
            if ( (*it)->marker() == boundMarker && (*it) != thisBoundary ){
                nextBoundary = (*it);
                break;
            }
        }
        if ( nextBoundary ){
            Node * nextNode = NULL; 
            if ( & nextBoundary->node( 0 ) != thisNode ){
                nextNode = & nextBoundary->node( 0 );
            } else {
                nextNode = & nextBoundary->node( 1 );
            }

            //** Check closed polygones here 
            if ( find( boundNodes.begin(), boundNodes.end(), nextNode ) == boundNodes.end() ) {
                boundNodes.push_front( nextNode );
            } else {
                break;
            }
        }
        thisBoundary = nextBoundary;
    }
    
    Mesh poly;
    std::vector < Node * > innerBound;
    for ( std::list < Node * >::iterator it = boundNodes.begin(); 
        it != boundNodes.end(); it ++ ){
        innerBound.push_back( poly.createNode( (*it)->pos() ) );
    }
        
    //** looking for polygon ends
    Node * upperLeftSurface = innerBound.front(); 
    Node * upperRightSurface = innerBound.back();
    if ( upperLeftSurface->pos()[ 0 ] > upperRightSurface->pos()[ 0 ] ){
        upperLeftSurface = innerBound.back(); 
        upperRightSurface = innerBound.front();
    }
    
    std::vector < Node * > outerBound;
    outerBound.push_back( upperLeftSurface );
    
    Node *n1 = poly.createNode( upperLeftSurface->pos() - RVector3( xBoundary, 0.0 ) );
    Node *n2 = poly.createNode( upperLeftSurface->pos() - RVector3( xBoundary, - mesh.ymin() + yBoundary ) );
    Node *n3 = poly.createNode( upperRightSurface->pos() - RVector3( -xBoundary, - mesh.ymin() + yBoundary ) );
    Node *n4 = poly.createNode( upperRightSurface->pos() - RVector3( -xBoundary, 0.0 ) );
    
    outerBound.push_back( upperLeftSurface );   
    RVector dx( increasingRange( 4.0, xBoundary, 20 ) );
    
    for ( uint i = 1; i < dx.size()-1; i ++ ) {
        outerBound.push_back( poly.createNode( upperLeftSurface->pos() 
                            + RVector3( -dx[ i ], 0 ) ) );
    }
    outerBound.push_back( n1 );   
    for ( uint i = 0; i < 9; i ++ ) {
        outerBound.push_back( poly.createNode( n1->pos() 
                            + ( n2->pos() - n1->pos() ) / 10 * ( i + 1 ) ) );
    }
    outerBound.push_back( n2 );   
    for ( uint i = 0; i < 9; i ++ ) {
        outerBound.push_back( poly.createNode( n2->pos() 
                            + ( n3->pos() - n2->pos() ) / 10 * ( i + 1 ) ) );
    }
    outerBound.push_back( n3 );   
    for ( uint i = 0; i < 9; i ++ ) {
        outerBound.push_back( poly.createNode( n3->pos() 
                            + ( n4->pos() - n3->pos() ) / 10 * ( i + 1 ) ) );
    }
    outerBound.push_back( n4 );   
    for ( uint i = dx.size()-2; i > 0; i -- ) {
        outerBound.push_back( poly.createNode( upperRightSurface->pos() 
                            + RVector3( dx[ i ], 0 ) ) );
    }
    outerBound.push_back( upperRightSurface );   
    
    for ( uint i = 0; i < outerBound.size() -1; i ++ ){
        poly.createEdge( *outerBound[ i ], *outerBound[ i + 1 ], -1 );
    }
    
    for ( uint i = 0; i < innerBound.size() -1; i ++ ){
        poly.createEdge( *innerBound[ i ], *innerBound[ i + 1 ], boundMarker );
    }
    
    Mesh boundaryMesh;
    TriangleWrapper( poly, boundaryMesh, "-YpzeAfaq" + str( 33 ) );
    
    std::vector < Boundary * > markedBoundaries2( boundaryMesh.findBoundaryByMarker( boundMarker ) );
    if ( markedBoundaries.size() == markedBoundaries2.size() ){
        for ( std::vector < Cell * >::const_iterator it = boundaryMesh.cells().begin(); 
                it != boundaryMesh.cells().end(); it ++ ){
            mesh.copyCell( *(*it) )->setMarker( cellMarker );
        }
    } else {
        return false;
        mesh.save("inMesh");
        boundaryMesh.save("boundaryMesh");
        throwError( 1, WHERE_AM_I + " Sorry!! Boundary mesh is not consistent to input mesh boundary. " 
                        + toStr( markedBoundaries.size() ) + " != " + toStr( markedBoundaries2.size() ) );
    }
    
    if ( save ){
        mesh.save("shortFOPinMesh");
        boundaryMesh.save("shortFOPboundaryMesh");
    }
    
    //** Fix boundary marker
    mesh.createNeighbourInfos( true );
    b = NULL;
    for ( std::vector < Boundary * >::const_iterator it = mesh.boundaries().begin();
        it != mesh.boundaries().end(); it ++ ){
        b = (*it);
        if ( ! b->leftCell() || ! b->rightCell() ) {
            if ( b->center().y() == mesh.ymax() ){
                b->setMarker( MARKER_BOUND_HOMOGEN_NEUMANN );
            } else {
                b->setMarker( MARKER_BOUND_MIXED );
            }
        }
    }
    return true;
}

} // namespace GIMLI
