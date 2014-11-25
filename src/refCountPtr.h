/***************************************************************************
 *   Copyright (C) 2014 by the resistivity.net development team            *
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

#ifndef _GIMLI_REFCOUNTPTR__H
#define _GIMLI_REFCOUNTPTR__H

#include "gimli.h"

namespace GIMLI{

// Base class for reference counted objects
// Scott Meyers: More Effective C++ Item 29 Source Code
class RCObject {                       
public:                                
    void addReference(){
        __M
        ++refCount_;
    }
  
    void removeReference(){
        __M
        if (--refCount_ == 0) delete this;
    }
    
    void markUnshareable(){
        shareable_ = false;
    }
    
    bool isShareable() const{
        return shareable_;
    }
    
    bool isShared() const {
        return refCount_ > 1;
    } 
 
protected:
    RCObject()
        : refCount_(0), shareable_(true) {}
   
    RCObject(const RCObject& rhs)
        : refCount_(0), shareable_(true) {}
   
    RCObject& operator=(const RCObject& rhs){
        return *this;
    }  
   
    virtual ~RCObject() = 0;

private:
    Index refCount_;
    bool shareable_;
};
 
 //*! Template class for smart pointers-to-T objects; 
 /*! Template class for smart pointers-to-T objects; 
  * must support the RCObject interface
  * */
template<class T > class RefCountIPtr {
public:      
     
    RefCountIPtr(T * realPtr = 0)
    : counter_(new CountHolder){ 
        counter_->pointee_ = realPtr;
        init_();
    }
   
    RefCountIPtr(const RefCountIPtr & rhs)
    : counter_(rhs.counter_){
        init_();
    }
   
    RefCountIPtr(RefCountIPtr & rhs)
    : counter_(rhs.counter_){
        THROW_TO_IMPL
    }
   
    ~RefCountIPtr(){
        counter_->removeReference();
    }
   
    RefCountIPtr & operator = (const RefCountIPtr & rhs){
        __M
        if (counter_ != rhs.counter_) {         
            counter_->removeReference();     
            counter_ = rhs.counter_;
            init_();
        }
        return *this;
    }
   
    RefCountIPtr & operator = (RefCountIPtr & rhs){
        __M
        if (counter_ != rhs.counter_) {         
            counter_->removeReference();     
            counter_ = rhs.counter_;
            init_();
        }
        THROW_TO_IMPL
        return *this;
    }
   
    T * operator->() {
        __M
        THROW_TO_IMPL
        return counter_->pointee_;
    }
    
    T * operator->() const{
        __M
        return counter_->pointee_;
    }
   
    T & operator*(){
        __M
        THROW_TO_IMPL
        return *counter_->pointee_;
    }
    
    T & operator*() const {
        __M
        return *counter_->pointee_;
    }
 
private:
    struct CountHolder : public RCObject {
        ~CountHolder() { delete pointee_; }
        T * pointee_;
    };
 
    void init_(){
        if (counter_->isShareable() == false) {
            T *oldValue = counter_->pointee_;
            counter_ = new CountHolder;
            counter_->pointee_ = oldValue ? new T(*oldValue) : 0;
        }     
        counter_->addReference();
    }
   
    CountHolder *counter_;
   
};
 
} // namespace GIMLI

#endif //_GIMLI_REFCOUNTPTR__H