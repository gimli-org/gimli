/***************************************************************************
 *   Copyright (C) 2009-2015 by the resistivity.net development team       *
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

#ifndef GIMLI_LONGOPTIONS__H
#define GIMLI_LONGOPTIONS__H

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>

#ifdef _MSC_VER
	#include <wingetopt.h>
#else
	#include <getopt.h>
#endif

#include "gimli.h"

namespace GIMLI{

class DLLEXPORT OptionBase{
public:
    virtual ~OptionBase(){}

    virtual void assign_(char * opt=NULL) = 0;

    virtual std::string typname() = 0;

    virtual std::string defautToString() = 0;

    inline const char key() const { return key_; }

    void assign(char * opt){
        if (hasArg_ == no_argument) assign_(NULL);
        if (hasArg_ == required_argument) assign_(opt);
    }

    void            * var_;
    char              key_;
    std::string       name_;
    std::string       help_;
    int               hasArg_;
};

inline std::ostringstream & operator << (std::ostringstream & s, const std::vector < std::string > & str){
    for (uint i = 0; i < str.size(); i ++) s << str[ i ] << " ";
    return s;
}

template < class T > class Option : public OptionBase{
public:
    Option(T & var){
        var_ = &var;
    }
    Option(T & var, const T & defaultVar){
        var_ = &var;
        *(T*)(var_) = defaultVar;
        defaultVar_ = defaultVar;
    }

    virtual ~Option(){}

    virtual void assign_(char * opt = NULL){ convert(*(T*)(var_), opt); }

    virtual std::string typname(){ return type(*(T*)(var_)); }

    T value() { return *(T*)(var_); }

    virtual std::string defautToString() {
        std::ostringstream streamOut;
        streamOut << defaultVar_;
        return streamOut.str();
    }

    T defaultVar_;
};

//! Simplified command line parser
/*! Simplified command line parser. The default option are -h [--help]. --version */
class DLLEXPORT OptionMap{
public:
    OptionMap(const std::string & description = "");

    virtual ~OptionMap();

    typedef std::map< char,         OptionBase * >    CharOptionMap;
    typedef std::map< std::string,  OptionBase * >    LongOptionMap;

    template < class T > void addLastArg(T & var, const std::string & lastArgString){
        lastArgString_ = lastArgString;
        lastOption_ = new Option< T >(var);
        lastOption_->hasArg_ = required_argument;
    }

    template < class T > void add(T & var, const std::string & key, const std::string & longOpt,
                                    const std::string & help){
        add(var, key, longOpt, help, T(var));
    }

    template < class T > void add(T & var, const std::string & key, const std::string & longOpt,
                                    const std::string & help, const T & defaultVar){

        if (key[ 0 ] != ':') allKeys_ += key;

        OptionBase * o = new Option< T >(var, defaultVar);
        o->name_ = longOpt;
        if ((*--key.end()) == ':' || (*--longOpt.end() ==':')) o->hasArg_ = required_argument;
        else o->hasArg_ = no_argument;

        std::string k( key.substr(0, key.rfind(':')));
        std::string lo(longOpt.substr(0, longOpt.rfind(':')));

        if (k.length() > 0) {
            o->key_ = k[ 0 ];
            cMap_.insert(std::pair< char, OptionBase * >(k[ 0 ], o));
        }

        lMap_.insert(std::pair< std::string, OptionBase * >(lo, o));
        options_.push_back(o);

        o->help_ = help;
    }

    void parse(int argc, char * argv[]);

    void buildLongOptions();

    inline void setDescription(const std::string & description){ descriptionString_ = description; }

    void printHelp(const std::string & main);

    OptionBase * findOption(const char & key) const {
//std::cout << "OptionBase * findOption(const char & key) const {" << std::endl;
        THROW_TO_IMPL
        std::list < OptionBase * >::iterator it;

        std::cout << std::equal_to< char >()(key, 'h') << std::endl;
        std::cout << std::equal_to< char >()(key, 'd') << std::endl;

        if (options_.size() > 0){
            std::cout << "tn: " << (*options_.begin())->typname() << std::endl;
            std::cout << std::mem_fun(&OptionBase::typname)(*options_.begin()) << std::endl;
            std::cout << "key: " << (*options_.begin())->key() << std::endl;
            std::cout << std::mem_fun(&OptionBase::key)(*options_.begin()) << std::endl;

            std::cout << "cmp: " << std::equal_to< char >()(key,
                                    std::mem_fun(&OptionBase::key)(*options_.begin()))
                      << std::endl;

            std::cout << "cmp: " << std::bind2nd(std::equal_to< char >(), key)('d') << std::endl;
            std::cout << "cmp: " << std::bind2nd(std::equal_to< char >(), key)('h') << std::endl;

            //** na schon fast, das geht leider noch nicht
//             std::cout << "cmp: " << std::bind2nd< std::mem_fun(&OptionBase::key) >(
//                     std::equal_to< char >(), key)(*options_.begin()) << std::endl;


// std::mem_fun(&OptionBase::key)))(*options_.begin())
//             << std::endl;

        }
//         it = std::find_if(options_.begin(), options_.end(),
//                         std::equal_to< char >()(key,
//                           std::mem_fun(&OptionBase::key)(*options_.begin())));
//                             std::equal_to< char > (std::mem_fun(&Option::key), key));

//std::cout << *it << std::endl;
/*
        if (it != options_.end()) {
            return *it;
        } else {
            return NULL;
        }*/
        return NULL;
    }

protected:
    struct option               * opts_;
    std::string                   allKeys_;
    std::string                   lastArgString_;
    std::string                   descriptionString_;

    OptionBase                  * lastOption_;
    std::list< OptionBase * >     options_;
    LongOptionMap                 lMap_;
    CharOptionMap                 cMap_;
    bool                          showVersion_;
};

} // namespace GIMLI
#endif //GIMLI_LONGOPTIONS__H
