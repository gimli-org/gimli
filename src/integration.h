/***************************************************************************
 *   Copyright (C) 2006-2017 by the GIMLi development team       *
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

#ifndef _GIMLI_INTEGRATION__H
#define _GIMLI_INTEGRATION__H

#include "gimli.h"

namespace GIMLI{

/*! */
class DLLEXPORT IntegrationRules : public Singleton< IntegrationRules >{
public:
    friend class Singleton< IntegrationRules >;
    
    /*! Return Gauss-Legendre quadrature point upto order <10. */
    inline const R3Vector & gauAbscissa(uint order) const { return gauAbscissa_[order]; }
    
    /*! Return Gauss-Legendre quadrature weights upto order <10. */
    inline const RVector & gauWeights(uint order) const { return gauWeights_[order]; }

    /*!
     * Generic quadrature positions for a triangle based on Gauss-Legendre quadrature
     H. T. RATHOD1*, K. V. NAGARAJA2, B. VENKATESUDU3 AND N. L. RAMESH4.
     Gauss Legendre quadrature over a triangle.
     J. Indian Inst. Sci., Sept.-Oct. 2004, 84, 183-188
    */
    inline const R3Vector & triGLAbscissa(uint order) const { return triGLAbscissa_[order]; }
    /*!
     * Generic quadrature weights for a triangle based on Gauss-Legendre quadrature
     H.T. RATHOD, K. V. NAGARAJA, B. VENKATESUDU AND N. L. RAMESH. 
     Gauss Legendre quadrature over a triangle.
     J. Indian Inst. Sci., Sept.-Oct. 2004, 84, 183-188
    */
    inline const RVector & triGLWeights(uint order) const { return triGLWeights_[order]; }

    inline const R3Vector & edgAbscissa(uint order) const { return edgAbscissa_[order]; }
    inline const RVector & edgWeights(uint order) const { return edgWeights_[order]; }

    inline const R3Vector & triAbscissa(uint order) const { return triAbscissa_[order]; }
    inline const RVector & triWeights(uint order) const { return triWeights_[order]; }

    inline const R3Vector & tetAbscissa(uint order) const { return tetAbscissa_[order]; }
    inline const RVector & tetWeights(uint order) const { return tetWeights_[order]; }

    inline const R3Vector & quaAbscissa(uint order) const { return quaAbscissa_[order]; }
    inline const RVector & quaWeights(uint order) const { return quaWeights_[order]; }

    inline const R3Vector & hexAbscissa(uint order) const { return hexAbscissa_[order]; }
    inline const RVector & hexWeights(uint order) const { return hexWeights_[order]; }

    inline const R3Vector & priAbscissa(uint order) const { return priAbscissa_[order]; }
    inline const RVector & priWeights(uint order) const { return priWeights_[order]; }
    /*! 
     * Return Gauss-Legendre quadrature positions for a given shape of the \ref MeshEntity upto order 10
     */
    const R3Vector & abscissa(const Shape & shape, uint order) const;
    /*!
     * Return Gauss-Legendre quadrature weights for a given shape of the \ref MeshEntity upto order 10
     */
    const RVector & weights(const Shape & shape, uint order) const;

    /*! Set whether triangle integration use Gauss Legendre polynomials (up to order 9) or 
     * native triangle coordinates (up to order 5). Default is GL == true. */
    inline void setTriUseGaussLegendre(bool use){ triUseGaussLegendre_ = use;}
    inline bool triUseGaussLegendre() const { return triUseGaussLegendre_; }
        
private:    
    /*! Default constructor */
    IntegrationRules();

    /*! Default destructor */
    ~IntegrationRules();

protected:
    
    void initGau_();
    void initTriGL_();

    void initEdg_();
    void initTri_();
    void initTet_();
    void initQua_();
    void initHex_();
    void initPri_();

    bool triUseGaussLegendre_;
    
    std::vector < R3Vector > gauAbscissa_;
    std::vector < RVector  > gauWeights_;

    /*!Quadrature positions for a triangle based on Gauss-Legendre quadrature */
    std::vector < R3Vector > triGLAbscissa_;
    /*!Quadrature weight for a triangle based on Gauss-Legendre quadrature */
    std::vector < RVector  > triGLWeights_;
    
    std::vector < R3Vector > edgAbscissa_;
    std::vector < RVector  > edgWeights_;
    
    std::vector < R3Vector > triAbscissa_;
    std::vector < RVector  > triWeights_;
    
    std::vector < R3Vector > tetAbscissa_;
    std::vector < RVector  > tetWeights_;
    
    std::vector < R3Vector > quaAbscissa_;
    std::vector < RVector  > quaWeights_;
    
    std::vector < R3Vector > hexAbscissa_;
    std::vector < RVector  > hexWeights_;
    
    std::vector < R3Vector > priAbscissa_;
    std::vector < RVector  > priWeights_;
    
};

} // namespace GIMLI{

#endif // _GIMLI_INTEGRATION__H
