/***************************************************************************
 *   Copyright (C) 2010 by the resistivity.net development team            *
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

#include <gimli.h>
#include <matrixTemplates.h>
#include <inversion.h>
#include <meshgenerators.h>
#include <string>
using namespace GIMLI;

#define vcout if ( verbose ) std::cout
#define dcout if ( debug ) std::cout
#define DEBUG if ( debug )

typedef KronMatrix< RMatrix, RMatrix > RKronMatrix;

class KronMatrixModelling : public ModellingBase {
public:
    KronMatrixModelling( RMatrix & A, RMatrix & B, bool verbose = false )
    : ModellingBase( verbose ), K_( A, B){
//        regionManager_->setParameterCount( K_.cols() );
        this->setMesh( createMesh1D( K_.cols() ) );
    }
    virtual ~KronMatrixModelling() { }
    
    RVector response( const RVector & model ) { return K_ * model; }

protected:
    RKronMatrix K_;
};

int main( int argc, char *argv [] ){

    RMatrix A, B;
    loadMatrixRow( A, "A.matrix" );
    loadMatrixRow( B, "B.matrix" );
    RVector x, y, z;
    load( x, "x.vector" );

    KronMatrixModelling f( A, B, true );

    y = f.response( x );
    std::cout << "y = " << y << std::endl;
//    f.regionManager().region( 0 )->setConstraintType ( 0 );

    RInversion inv( y, f, true, false );
    inv.setModel( RVector( x.size(), 1.0 ) );
    inv.setAbsoluteError ( 0.1 );
    z = inv.run();
    std::cout << "z = " << z << std::endl;

    return EXIT_SUCCESS;
}
