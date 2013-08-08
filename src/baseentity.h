/***************************************************************************
 *   Copyright (C) 2007-2011 by the resistivity.net development team       *
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

#ifndef _GIMLI_BASEENTITY__H
#define _GIMLI_BASEENTITY__H 

#include "gimli.h"

namespace GIMLI{

  //! Base Entity
  /*! The is a base entity which holds some basic informations. An identification 
  * number for the object, a validation status flag and a runtime type information (RTTI). 
  * This information will be inherited in many objects i.e. the mesh and data entitys. */
class DLLEXPORT BaseEntity{
    public:
    BaseEntity() : id_( -1 ), valid_( false ), marker_( 0 ) { }

    BaseEntity( const BaseEntity & ent ) 
        : id_( ent.id() ), valid_( ent.valid() ), marker_( ent.marker() ){ 
    }

    BaseEntity & operator = ( const BaseEntity & ent ){
        if ( this != &ent ){
            id_ = ent.id();
            valid_ = ent.valid();
            marker_ = ent.marker();
        } return * this;
    }
  
    virtual ~BaseEntity( ){}

    /*! Return entity rtti value. */
    inline virtual uint rtti() const { return MESH_BASEENTITY_RTTI; }

    /*! Return entity valid status. */
    inline virtual bool valid() const { return valid_; }

    /*! Set entity valid status. */
    inline virtual void setValid( bool valid ) { valid_ = valid ; }

    /*! Returns the entity id. */
    inline int id() const { return id_; }

    /*! Set the entity id. */
    inline void setId( int id ) { id_ = id ; }
    
    inline void setMarker( int marker ) { marker_ = marker; }
    
    inline int marker() const { return marker_; }
    
protected:
    int id_;
    
    bool valid_;
    
    int marker_;
}; // class BaseEntity;

}  // namespace GIMLI{


#endif // _GIMLI_BASEENTITY__H
