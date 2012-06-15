/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
 *   Thomas Günther thomas@resistivity.net                                 *
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

#ifndef _GIMLI_TRANS__H
#define _GIMLI_TRANS__H

#include "gimli.h"
#include "vector.h"
#include "numericbase.h"

#define TRANSTOL 1e-8

namespace GIMLI{

/*! Basis Transform vector (unity), Identity transformation,
that means it do nothing. Just for prototyping f(a). */
template< class Vec > class Trans {
public:
//         
    /*! Default constructur. */
    Trans( ) { }

    /*! Default destructor.*/
    virtual ~Trans( ) { }

    /*! Return a. */
    virtual Vec trans( const Vec & a ) const { return a; }

    /*! */
    virtual Vec invTrans( const Vec & a ) const { return a; }

    /*! */
    virtual Vec deriv( const Vec & a ) const { return Vec( a.size(), 1.0 ); }

    /*! Update parameter by df: invTrans( f( a ) + df ). \n
    intrinsic function that have never to be overloaded. */
    Vec update( const Vec & a, const Vec & b ) const {
        return ( invTrans( trans( a ) + b ) );
    }

    /*! Error of f(a) calculated by a and relative error da/a: df = | da * df/da | \n
    intrinsic function that could be overloaded */
    Vec error( const Vec & a, const Vec & daBya ) const {
        return abs( Vec( a * daBya * deriv( a ) ) );
    }

    /*! Alternative Version - brute force: df = | f( a + da ) - f( a) | \n
    intrinsic functions that have never to be overloaded */
    Vec error_brute( const Vec & a, const Vec & daBya ) const {
        return abs( trans( a * ( daBya + 1.0 ) ) - trans ( a ) );
    }
};

/*! Base class for non-invertible transformations, e.g. transMult and transPlus \n
    making inverse transform using Newton's method */
template< class Vec > class TransNewton : public Trans < Vec >{
public:
    TransNewton( const int maxiter = 10 ) : maxiter_( maxiter ) { }
    virtual ~TransNewton( ){ }

    virtual Vec trans( const Vec & a ) const { }

    virtual Vec deriv( const Vec & a ) const { }

    virtual Vec invTrans( const Vec & a ) const {
        Vec p( Vec( a.size(), 0.01 ) ); // guess vector
        Vec dm( trans( p ) - a ); // function to be nulled
        //std::cout << "Newton iter norm:";
        for( size_t i = 0; i < 10; i++ ) {
            p = p - dm / deriv( p ); // Newton method
            dm = trans( p ) - a;
            //std::cout << " " << norml2( dm );
        }
        //std::cout << std::endl;
        return p;
    }
protected:
    int maxiter_;
};

/*! Fundamental transformation operations: linear, inverse, log, power, exp: */
/*! Linear multiplication (e.g. geometric factor) and addition f(a) = a * factor + offset. */
template< class Vec > class TransLinear : public Trans < Vec > {
public:
    TransLinear( const Vec & factor, const Vec & offset ) : factor_( factor ), offset_( offset ){ }
    TransLinear( const Vec & factor, double offset = 0.0 ) : factor_( factor ), offset_( Vec( factor.size(), offset ) ) {}
    virtual ~TransLinear( ) { }

    virtual Vec trans( const Vec & a ) const {
//        std::cout << a.size() << " " << factor_.size() << " " << offset_.size()<< std::endl;
        return a * factor_ + offset_;
    }
    virtual Vec invTrans( const Vec & a ) const { return ( a - offset_ ) / factor_; }
    virtual Vec deriv( const Vec & a ) const { return factor_; }

protected:
    Vec factor_;
    Vec offset_;
};

/*! power transformation f(a) = pow( a / a_0 , n ). */
template< class Vec > class TransPower : public Trans < Vec > {
public:
    TransPower( const double npower = -1.0, const double a0 = 1.0 ) : npower_( npower ), a0_( a0 ) { }
    virtual ~TransPower( ) { }

    virtual Vec trans( const Vec & a ) const { return pow( a / a0_ , npower_ ); }
    virtual Vec invTrans( const Vec & f ) const { return pow( f , 1.0 / npower_ ) * a0_; }
    virtual Vec deriv( const Vec & a ) const { return pow( a / a0_ , npower_ - 1.0 ) * npower_ / a0_; }
protected:
    double npower_;
    double a0_;
};

/*! Exponential transformation f(a) = f_0 exp( - a / a_0 ). */
template< class Vec > class TransExp : public Trans < Vec > {
public:
    TransExp( const double a0 = 1.0, const double f0 = 1.0 ) : a0_( a0 ), f0_( f0 ) { }
    virtual ~TransExp( ) { }

    virtual Vec trans( const Vec & a ) const { return exp( a / a0_ * (-1.0) ) * f0_; }
    virtual Vec invTrans( const Vec & f ) const { return log( f / f0_ ) * a0_ * (-1.0); }
    virtual Vec deriv( const Vec & a ) const { return trans( a ) / a0_ * (-1.0); }
protected:
    double a0_;
    double f0_;
};

/*! Inverse of the parameter (e.g. velocity/slowness) f(a) = 1 / a. */
template< class Vec > class TransInv : public Trans < Vec > {
public:
    TransInv( ) { }
    virtual ~TransInv( ) { }

    virtual Vec trans( const Vec & a ) const { return 1.0 / a; }
    virtual Vec invTrans( const Vec & f ) const { return 1.0 / f; }
    virtual Vec deriv( const Vec & a ) const { return -1.0 / a / a; }

};

/*! Logarithm of the parameter with the natural bound 0 or a self-defined. */
template< class Vec > class TransLog : public Trans < Vec > {
public:
    TransLog( double lowerbound = 0.0 ) : lowerbound_( lowerbound ) { }
    virtual ~TransLog( ) { }
    
    virtual Vec trans( const Vec & a ) const {
        double lb1 = lowerbound_ * ( 1.0 + TRANSTOL );
        if ( min( a ) < lb1 ) {
            std::cerr << WHERE_AM_I << " Warning! " << min( a )
                      << " <=" << lowerbound_ << " lowerbound" << std::endl;
            Vec tmp( a );
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = max( a[ i ], lb1  );
            }
            return log( tmp - lowerbound_ );
        }
        return log( a - lowerbound_ );
    }

    virtual Vec invTrans( const Vec & f ) const { return exp( f ) + lowerbound_; }

    virtual Vec deriv( const Vec & a ) const { 
        double lb1 = lowerbound_ * ( 1.0 + TRANSTOL );
        if ( min( a ) < lb1 ) {
            std::cerr << WHERE_AM_I << " Warning! " << min( a )
                      << " <=" << lowerbound_ << " lowerbound" << std::endl;
            Vec tmp( a );
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = max( a[ i ], lb1  );
            }
            return 1.0 / ( tmp - lowerbound_ ); 
        }
        return 1.0 / ( a - lowerbound_ ); 
    }

    //** suggested by Friedel(2003), but deactivated since inherited by transLogLU 
//    virtual Vec error( const Vec & a, const Vec & daBya ) const { return log( 1.0 + daBya ); }

    inline void setLowerBound( double lb ) { lowerbound_ = lb; }

    inline double lowerBound() const { return lowerbound_; }

protected:
    double lowerbound_;
};


/*! Range constraint transform functions: */
/*! Logarithmic barrier with upper and lower bound f(a) = log( a - a_u ) - log( a_u - a ). */
template< class Vec > class TransLogLU : public TransLog < Vec > {
public:
    TransLogLU( double lowerbound = 0.0, double upperbound = 0.0 )
        : TransLog< Vec >( lowerbound ), upperbound_( upperbound ) { }
    virtual ~TransLogLU( ) { }

    Vec rangify( const Vec & a ) const {
        Vec tmp( a );
        double lb1 = this->lowerBound() * ( 1.0 + TRANSTOL );
        if ( min( a ) < lb1 ){
            std::cerr << WHERE_AM_I << " Warning! " << min( a )
                      << " <=" << this->lowerBound() << " lower bound" << std::endl;
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = max( a[ i ], lb1 );
            }
        }

        double ub1 = upperbound_ * ( 1.0 - TRANSTOL );
        if ( max( a ) > ub1  ){
            std::cerr << WHERE_AM_I << " Warning! " << max( a ) << " > "
                      << upperbound_ << " upper bound" << std::endl;
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = min( tmp[ i ], ub1 );
            }
        }
        return tmp;
    }
    
    virtual Vec trans( const Vec & a ) const {
        if ( std::fabs( upperbound_ ) < TOLERANCE ) return TransLog< Vec >::trans( a );

        Vec tmp = rangify( a );
        return ( log( tmp - this->lowerBound() ) - log( upperbound_ - tmp ) );
    }

    virtual Vec invTrans( const Vec & a ) const {
        if ( std::fabs( upperbound_ ) < TOLERANCE ) return TransLog< Vec >::invTrans( a );

        //return ( ( exp( a ) * upperbound_ + this->lowerBound() ) / ( exp( a ) + this->lowerBound() ) ); //** keine Ahnung
        return ( ( exp( a ) * upperbound_ + this->lowerBound() ) / ( exp( a ) + 1.0 ) );
    }

    virtual Vec deriv( const Vec & a ) const {
        if ( std::fabs( upperbound_ ) < TOLERANCE ) return TransLog< Vec >::deriv( a );

        Vec tmp = rangify( a );
        return ( 1.0 / ( tmp - this->lowerBound() ) + 1.0 / ( upperbound_ - tmp ) );
    }

protected:
  double upperbound_;
};

/*! Cotangens barrier method, e.g. for water content (NMR) */
template< class Vec > class TransCotLU : public Trans < Vec > {
public:

    TransCotLU( double lowerbound = 0.0, double upperbound = 0.0 )
        : lowerbound_( lowerbound ), upperbound_( upperbound ) { }
        
    virtual ~TransCotLU( ) { }

    virtual Vec trans( const Vec & a ) const {
        Vec tmp( a );
        double fak = 1.00001;
        if ( min( a ) <= lowerbound_  ){
            std::cerr << WHERE_AM_I << " Warning! " << min( a ) << " < " << lowerbound_ << " = lowerbound" << std::endl;
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = max( a[ i ], lowerbound_ * fak );
            }
        }
        if ( max( a ) >= upperbound_  ){
            std::cerr << WHERE_AM_I << " Warning! " << max( a ) << " > " << upperbound_ << " = upperbound" << std::endl;
            for ( uint i = 0; i < a.size(); i ++ ){
                tmp[ i ] = min( a[ i ], upperbound_ / fak );
            }
        }
        return cot( ( tmp - lowerbound_ ) / ( upperbound_ - lowerbound_ ) * PI ) * (-1.0);
    }

    virtual Vec invTrans( const Vec & a ) const {
//    return acot( a * (-1.0) ) * ( upperbound_ - lowerbound_ ) / PI + lowerbound_;
        return atan( a ) * ( upperbound_ - lowerbound_ ) / PI + ( lowerbound_ + upperbound_ ) / 2;
    }

    virtual Vec deriv( const Vec & a ) const {
        return ( ( trans( a ) * trans( a ) + 1.0 ) * PI / ( upperbound_ - lowerbound_ ) );
    }

protected:
    double lowerbound_;
    double upperbound_;
};

/*! Different combined transformations: */
/*! Two nested TF's, e.g. log( G * resistance ) or LogLU( CRIM( slowness ) ) */
template< class Vec > class TransNest : public Trans < Vec >{
public:
    TransNest( Trans< Vec > & T1, Trans< Vec > & T2 ) : T1_( &T1 ), T2_( &T2 ){
    }
    virtual ~TransNest( ){ }

    virtual Vec trans( const Vec & a ) const {
        return T1_->trans( T2_->trans( a ) );
    }

    virtual Vec invTrans( const Vec & a ) const {
        return T2_->invTrans( T1_->invTrans( a ) );
    }

    virtual Vec deriv( const Vec & a ) const {
        return T1_->deriv( T2_->trans( a ) ) * T2_->deriv( a );
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Two added TF's f(a) = f1(a) + f2(a) */
template< class Vec > class TransAdd : public TransNewton < Vec >{
public:
    TransAdd( Trans< Vec > & T1, Trans< Vec > & T2 ) : T1_( &T1 ), T2_( &T2 ){ }
    virtual ~TransAdd( ){ }

    virtual Vec trans( const Vec & a ) const {
        return T1_->trans( a ) + T2_->trans( a );
    }

    virtual Vec deriv( const Vec & a ) const {
        return T1_->deriv( a ) + T2_->deriv( a );
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Two multiplied TF's f(a) = f1(a) * f2(a) (not necessarily stable). */
template< class Vec > class TransMult : public TransNewton < Vec >{
public:
    TransMult( Trans< Vec > & T1, Trans< Vec > & T2 ) : T1_( &T1 ), T2_( &T2 ){ }
    virtual ~TransMult( ){ }

    virtual Vec trans( const Vec & a ) const {
        return T1_->trans( a ) * T2_->trans( a );
    }

    virtual Vec deriv( const Vec & a ) const {
        return T1_->deriv( a) * T2_->trans( a ) + T2_->deriv( a ) * T1_->trans( a );
    }

protected:
    Trans< Vec > * T1_;
    Trans< Vec > * T2_;
};

/*! Petrophysical transform functions (specializations&combinations of other) */
/*! Complex refractive index model (CRIM) based on transLinear */
template< class Vec > class TransCRIM : public TransLinear < Vec > {
public:
    /*! petrophysical constructor using porosity and permittivity of matrix and water as input */
    TransCRIM( const Vec & phi, double ematrix = 5.0, double ewater = 81.0 )
        : TransLinear< Vec >( 1.0 / phi / 2.99e8 / ( std::sqrt( ewater ) - 1.0 ),
          ( phi + ( phi - 1.0 ) * std::sqrt( ematrix ) ) / ( phi * ( std::sqrt( ewater ) - 1.0 ) ) ){
    }
    virtual ~TransCRIM( ) { }
};

/*! Archie's law transforming resistivity into porosity. */
template< class Vec > class TransArchie : public TransPower < Vec > {
public:
    /*! petrophysical constructor using fluid resistivity as input */
    TransArchie( double resfluid, double mexp = 2.0 )
        :TransPower< Vec >( resfluid, 1 / mexp ){
    }

    virtual ~TransArchie( ) { }
};

/*! obsolete transform functions, e.g. for test reasons */
/*! transQuadrad, test class for checking transNewton */
template< class Vec > class TransQuadrat : public TransNewton < Vec >{
public:
    TransQuadrat( const int maxiter = 10 ) : TransNewton< Vec >( maxiter ){ }
    virtual ~TransQuadrat( ){ }
    virtual Vec trans( const Vec & a ) const { return a * a; }
    virtual Vec deriv( const Vec & a ) const { return a * 2.0; }
};

/*! Logarithm of parameter times factor (e.g. log( G * resistance )), better by transNest(transLog,transLinea) */
template< class Vec > class TransLogMult : public TransLog < Vec > {
public:
    TransLogMult( const Vec & factor, double lowerbound = 0.0 )
        : TransLog< Vec >( lowerbound ), factor_( factor ) { }
    virtual ~TransLogMult( ) { }

    virtual Vec trans( const Vec & a ) const { return TransLog< Vec >::trans( a * factor_ ); }
    virtual Vec invTrans( const Vec & a ) const { return TransLog< Vec >::invTrans( a ) / factor_; }
    virtual Vec deriv( const Vec & a ) const { return TransLog< Vec >::deriv( a * factor_ ) * factor_; }

protected:
    Vec factor_;
};

/*! Tangens barrier method, replaced by TransCotLU, same results */
template< class Vec > class TransTanLU : public Trans < Vec > {
public:
    TransTanLU( double lowerbound = 0.0, double upperbound = 0.0 )
        : lowerbound_( lowerbound ), upperbound_( upperbound ) {
    }
    virtual ~TransTanLU( ) { }

    virtual Vec trans( const Vec & a ) const {
        Vec tmp( a );
        double tol = 1e-5;
        if ( a.min() < lowerbound_  ){
            std::cerr << WHERE_AM_I << " Warning! " << a.min() << " < " << lowerbound_ << " lowerbound" << std::endl;
            for ( int i = 0; i < a.size(); i ++ ){
                tmp[ i ] = max( a[ i ], lowerbound_ ) + tol;
            }
        }
        if ( a.max() > upperbound_  ){
            std::cerr << WHERE_AM_I << " Warning! " << a.max() << " > " << upperbound_ << " lowerbound" << std::endl;
            Vec tmp( a );
            for ( int i = 0; i < a.size(); i ++ ){
                tmp[ i ] = min( a[ i ], upperbound_ ) - tol;
            }
        }
        return tan( ( ( a - lowerbound_ ) / ( upperbound_ - lowerbound_ ) - 0.5 ) * PI );
    }

    virtual Vec invTrans( const Vec & a ) const {
        return atan( a ) / PI * ( upperbound_ - lowerbound_ ) + ( lowerbound_ + upperbound_ ) / 2;
    }

    virtual Vec deriv( const Vec & a ) const {
        return ( pow( trans( a ), 2.0 ) + 1.0 );
    }

protected:
    double lowerbound_;
    double upperbound_;
};

/*! Logarithm barrier of parameter times factor, better by transNest(transLogLU,transLinear) */
template< class Vec > class TransLogLUMult : public TransLogLU < Vec > {
public:
    TransLogLUMult( const Vec & factor, double lowerbound = 0.0, double upperbound = 0.0 )
        : TransLogLU< Vec >( lowerbound, upperbound ), factor_( factor ){ }
    virtual ~TransLogLUMult( ) { }

    virtual Vec trans( const Vec & a ) const {
        return TransLogLU< Vec >::trans( a * factor_ );
    }
    virtual Vec invTrans( const Vec & a ) const {
        return TransLogLU< Vec >::invTrans( a ) / factor_ ;
    }
    virtual Vec deriv( const Vec & a ) const {
        return TransLogLU< Vec >::deriv( a * factor_ ) * factor_;
    }

protected:
    Vec factor_;
};

/*! Auxillary cumulative transformation function using a vector of transformations*/
//! Very Slow. Refactor it!!
template< class Vec > class CumulativeTrans : public Trans <Vec>{
public:
    CumulativeTrans( ) { }

    virtual ~CumulativeTrans( ) { }

    virtual Vec trans( const Vec & a ) const {
//         Vec tmp;
//         //std::cout << "CumulativeTrans::Trans" << std::endl;
//         for ( uint i = 0; i < transVec_.size(); i ++ ){
//             Vec ap( a, bounds_[ i ].first, bounds_[ i ].second );
//             //std::cout << i << " " << bounds_[ i ].first << ": " << bounds_[ i ].second << std::endl;
// //              std::cout << tmp << std::endl;
//             tmp = cat( tmp, transVec_[ i ]->trans( ap ) );
//           //  std::cout << tmp << std::endl;
//         }
//         return tmp;
        Vec tmp( a.size() );
        for ( uint i = 0; i < transVec_.size(); i ++ ){
            tmp.setVal( transVec_[ i ]->trans( a(bounds_[ i ].first, bounds_[ i ].second) ),
                        bounds_[ i ].first, bounds_[ i ].second );
        }
        return tmp;
    }

    virtual Vec invTrans( const Vec & a ) const {
        Vec tmp( a.size() );
        for ( uint i = 0; i < transVec_.size(); i ++ ){
            tmp.setVal( transVec_[ i ]->invTrans( a(bounds_[ i ].first, bounds_[ i ].second) ),
                        bounds_[ i ].first, bounds_[ i ].second );
        }
        return tmp;
//         Vec tmp;
//         for ( uint i = 0; i < transVec_.size(); i ++ ){
//             Vec ap( a, bounds_[ i ].first, bounds_[ i ].second );
//             tmp = cat( tmp, transVec_[ i ]->invTrans( ap ) );
//         }
//         return tmp;
    }

    virtual Vec deriv( const Vec & a ) const {
//         Vec tmp;
//         for ( uint i = 0; i < transVec_.size(); i ++ ){
//             Vec ap( a, bounds_[ i ].first, bounds_[ i ].second );
//             tmp = cat( tmp, transVec_[ i ]->deriv( ap ) );
//         }
        Vec tmp( a.size() );
        for ( uint i = 0; i < transVec_.size(); i ++ ){
            tmp.setVal( transVec_[ i ]->deriv( a(bounds_[ i ].first, bounds_[ i ].second) ),
                        bounds_[ i ].first, bounds_[ i ].second );
        }
        return tmp;
    }
    void push_back( Trans< Vec > & trans , uint len ) {
        uint oldlen = 0;
        if ( ! bounds_.empty() ) oldlen = bounds_.back().second;
        transVec_.push_back( &trans );
        bounds_.push_back( std::pair< uint, uint >( oldlen, oldlen + len ) );
    }

    std::vector < Trans< Vec > * > transVec_;
    std::vector < std::pair< uint, uint > > bounds_;
};

typedef Trans < RVector > RTrans;
typedef TransLinear < RVector > RTransLinear;
typedef TransPower < RVector > RTransPower;
typedef TransLog < RVector > RTransLog;
typedef TransLogLU < RVector > RTransLogLU;
typedef TransCotLU < RVector > RTransCotLU;
typedef CumulativeTrans< RVector > RCumulativeTrans;
} // namespace GIMLI


#endif // _GIMLI_TRANS__H

