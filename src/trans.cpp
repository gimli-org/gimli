/******************************************************************************
 *   Copyright (C) 2006-2019 by the GIMLi development team                    *
 *   Thomas Günther thomas@resistivity.net                                    *
 *   Carsten Rücker carsten@resistivity.net                                   *
 *                                                                            *
 *   Licensed under the Apache License, Version 2.0 (the "License");          *
 *   you may not use this file except in compliance with the License.         *
 *   You may obtain a copy of the License at                                  *
 *                                                                            *
 *       http://www.apache.org/licenses/LICENSE-2.0                           *
 *                                                                            *
 *   Unless required by applicable law or agreed to in writing, software      *
 *   distributed under the License is distributed on an "AS IS" BASIS,        *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 *   See the License for the specific language governing permissions and      *
 *   limitations under the License.                                           *
 *                                                                            *
 ******************************************************************************/

#include "trans.h"

namespace GIMLI{

template <> DLLEXPORT RVector TransLogLU<RVector>::rangify(const RVector & a) const {
    RVector v(a);
    double minA = min(v);
    double maxA = max(v);
    double lb1 = this->lowerBound() * (1.0 + TRANSTOL);
    if (minA < lb1){
        // log(Verbose, "TransLogLU cap min: ",
        //                 minA, "< lower bound(", this->lowerBound(),")");
        capMin(v, lb1);
    }

    double ub1 = upperbound_ * (1.0 - TRANSTOL);
    if (maxA > ub1 ){
        // log(Verbose, "TransLogLU cap max: ",
        //                 maxA, "> upper bound(", upperbound_,")");
        capMax(v, ub1);
    }
    return v;
}

template <> RVector TransLogLU<RVector>::trans(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::trans(a);

    RVector tmp = rangify(a);
    return (log(tmp - this->lowerBound()) - log(upperbound_ - tmp));
}

template <> RVector TransLogLU<RVector>::invTrans(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::invTrans(a);
    RVector expm(a);
    capMax(expm, 50.0);
    expm = exp(expm);
    return (expm * this->upperBound() + this->lowerBound()) / (expm + 1.0);
}

template <> RVector TransLogLU<RVector>::deriv(const RVector & a) const {
    if (std::fabs(upperbound_) < TOLERANCE) return TransLog< RVector >::deriv(a);
    RVector tmp = rangify(a);
    return (1.0 / (tmp - this->lowerBound()) + 1.0 / (this->upperBound() - tmp));
}

} // namespace GIMLI
