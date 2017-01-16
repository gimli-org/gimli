/***************************************************************************
 *   Copyright (C) 2012-2014 by the GIMLi development team       *
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

#ifndef _GIMLI_GRAVIMETRY__H
#define _GIMLI_GRAVIMETRY__H

#include "gimli.h"
#include "modellingbase.h"

namespace GIMLI {

//! Modelling class for gravimetry calculation using polygon integration
/*! Modelling class for gravimetry calculation using polygon integration */
class DLLEXPORT GravimetryModelling : public ModellingBase {
public:
    GravimetryModelling( Mesh & mesh, DataContainer & dataContainer, bool verbose = false );

    virtual ~GravimetryModelling() { }

    RVector createDefaultStartModel( );

    /*! Interface. Calculate response */
    virtual RVector response( const RVector & slowness );

    /*! Interface. */
    virtual void createJacobian( const RVector & slowness );

    /*! Interface. */
    virtual void initJacobian( );

protected:
};



/*! Only for a small TPOC for 2d gravimetry after WonBevis1987
 Do not use until u know what u do. */
DLLEXPORT double lineIntegraldGdz( const RVector3 & p1, const RVector3 & p2 );

/*! Do not use until u know what u do. */
DLLEXPORT RVector calcGBounds( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model );

/*! Do not use until u know what u do. */
DLLEXPORT RVector calcGCells( const std::vector< RVector3 > & pos, const Mesh & mesh, const RVector & model, uint nInt = 0 );

} //namespace GIMLI

#endif
