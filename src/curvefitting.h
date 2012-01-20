/***************************************************************************
 *   Copyright (C) 2009-2011 by the resistivity.net development team       *
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

#ifndef _GIMLI_CURVEFITTING__H
#define _GIMLI_CURVEFITTING__H

#include "gimli.h"
#include "modellingbase.h"

namespace GIMLI{

/*! Abstract function class, i.e., functor */
template < class ArgType , class ValueType > class DLLEXPORT Function {
public:
    Function( ) { }
    
    virtual ~Function( ) { }
    
    Function( const Function & funct ){ copy_( funct ); }
    
    Function & operator = ( const Function & funct ){
        if ( this != & funct ){
            copy_( funct );
        } return * this;
    }
    
    inline ValueType operator() ( const ArgType & arg ) const { return this->getValue( arg ); }
    
    inline Vector < ValueType > operator() ( const Vector < ArgType > & arg ) const { return this->getValue( arg ); }
    
    virtual ValueType getValue( const ArgType & arg ) const = 0;
    
    virtual Vector < ValueType > getValue( const Vector < ArgType > & arg ) const = 0;
    
protected:
    virtual void copy_( const Function & funct ) { };
};

class DLLEXPORT HarmonicFunction : public Function< double, double >{
public:
    HarmonicFunction( const RVector & coeff, double xmin, double xmax );
    
    virtual ~HarmonicFunction( );

    virtual double getValue( const double & arg ) const;
    
    virtual RVector getValue( const RVector & arg ) const;
    
    void setCoefficients( const RVector & coeff );
    
    inline const RVector & coefficients() const { return coeff_; }
    
    inline void setXMin( double xmin ){ xMin_ = xmin; }
    
    inline double xMin() const { return xMin_; }
    
    inline void setXMax( double xmax ){ xMax_ = xmax; }

    inline double xMax() const { return xMax_; }
    
protected:
    virtual void copy_( const HarmonicFunction & funct );

    RVector coeff_;
    size_t nHarmonic_;
    double xMin_;
    double xMax_;
};

class DLLEXPORT HarmonicModelling : public ModellingBase {
public:
    /*! constructor, nh: number of coefficients, xvec: abscissa, */
    HarmonicModelling( size_t nh, const RVector & tvec, bool verbose = false  );
    
    virtual ~HarmonicModelling( ){ }
    
    /*! the main thing - the forward operator: return f(x) */
    virtual RVector response( const RVector & par );

    /*! an additional forward operator for another time basis */
    virtual RVector response( const RVector & par, const RVector tvec );

    /*! optional: generation of jacobian matrix, uncomment for default behaviour (brute force) */
    virtual void createJacobian( const RVector & model );

    /*! Define the startmodel */
    inline virtual RVector startModel( ){ return RVector( np_, 0.0 ); }
    
protected:
    
    RVector t_; //! abscissa vector x
    RMatrix A_; //! function matrix
    double tMin_;
    double tMax_;
    size_t nh_, nt_, np_;
};


} // namespace GIMLI{

#endif // _GIMLI_CURVEFITTING__H
