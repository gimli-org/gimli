/***************************************************************************
 *   Copyright (C) 2009-2011     by the resistivity.net development team   *
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

#ifndef _GIMLI_INVERSION_FRAMEWORK__H
#define _GIMLI_INVERSION_FRAMEWORK__H

#include "gimli.h"

// This may be the future, just collect ideas and concepts here.

//! Basic framework to perform an inversion
/*!*/
class InversionFramework {
public:
    /*!*/
    InversionFramework( Inversion & inversion ) : inversion_( & inversion ){
    }

    /*! Default destructor */
    virtual InversionFramework(){}

protected:
    Inversion * inversion_;
}

//! Framework to perform a time lapse inversion
/*!*/
class InversionTimelapse : public InversionFramework  {
public:
    /*! Default destructor */
    virtual InversionFramework(){}
protected:
};

//! Framework to perform a roll along inversion
/*!*/
class InversionRollalong : public InversionFramework {
public:
    /*! Default destructor */
    virtual InversionFramework(){}
protected:    
}
    
//! Framework to perform a roll along inversion, along the time dimension
/*!*/
class InversionRollalongInTime : public InversionRollalong {
public:
    /*! Default destructor */
    virtual InversionFramework(){}
protected:    
}

//! Framework to perform a roll along inversion, along the space dimension
/*!*/
class InversionRollalongInSpace : public InversionRollalong {
public:
    /*! Default destructor */
    virtual InversionFramework(){}
protected:    
};

