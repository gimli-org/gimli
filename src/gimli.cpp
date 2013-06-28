/***************************************************************************
 *   Copyright (C) 2006-2007 by the resistivity.net development team       *
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

#include "gimli.h"

#include <cerrno>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace GIMLI{

//extern bool __SAVE_PYTHON_GIL__ = false;
//extern bool __GIMLI_DEBUG__ = false;

std::string authors(){
  std::string a( (std::string)("bugs and suggestions to: ") + PACKAGE_AUTHORS );
  return a;
}

int openFile(const std::string & fname, std::fstream * file, std::ios_base::openmode farg, bool terminate){
    file->open(fname.c_str(), farg);
    if (!*file){
        if (terminate) {
            throwError( EXIT_OPEN_FILE, WHERE_AM_I + " '" + fname + "': " +strerror(errno) + str(errno));
        }
        return 0;
    }
    return 1;
}

bool fileExist( const std::string & filename ){
    bool result = false;
    std::ifstream file; file.open( filename.c_str() );
    if ( file ) {
        result = true;
        file.close();
    }
    return result;
}

uint fileLength( std::fstream & file ){
    uint oldPos = file.tellg();
    file.seekg( 0, std::ios::end );
    uint length = file.tellg();
    file.seekg( oldPos );
    return length;
}

uint countColumnsInFile( const std::string & fname){
    uint columnsCount = 0;
    return countColumnsInFile( fname, columnsCount );
}

uint countColumnsInFile( const std::string & fname, uint & columnsCount ){
    columnsCount = 0;
    std::fstream file; if ( !openInFile( fname, & file, false ) ) { return 0; }

    std::vector < std::string > subStrings;
    std::string str, tmp;

    while ( !file.eof() ){
        getline( file, str );
        if ( str.find( '#' ) != std::string::npos || str.empty() ) {
            columnsCount++;
        } else {
            file.close();
            return getSubstrings( str ).size();
        }
    }

    file.close();
    return 0;
}

uint countRowsInFile( const std::string & fname ){
    std::fstream file; openInFile( fname, & file );
    CERR_TO_IMPL;
    file.close();
    return 0;
}

std::vector < std::string > getRowSubstrings( std::fstream & file, char comment ){
  std::vector < std::string > subStrings;
  std::string str, tmp;
  getline( file, str );

  std::istringstream is( str.substr( 0, str.find( comment ) ) );
  while ( is >> tmp ) subStrings.push_back( tmp );
  return subStrings;
}

std::vector < std::string > getNonEmptyRow( std::fstream & file, char comment ){
  std::vector < std::string > row;
  while( (row = getRowSubstrings( file, comment ) ).empty() && !file.eof() );
  return row;
}

std::vector < std::string > getCommentLine( std::fstream & file, char comment ){
    std::vector< std::string > row;
    std::string str;
    getline( file, str );
    row = getSubstrings( str.substr( str.find( comment ), std::string::npos ) );
    return row;
}

std::vector < std::string > getSubstrings( const std::string & str ){
  std::vector < std::string > subStrings;
  std::istringstream is( str );
  std::string tmp;
  while ( is >> tmp ) subStrings.push_back( tmp );
  return subStrings;
}

std::vector < std::string > split( const std::string & str, char delimiter ){
    std::vector < std::string > subStrings;
    size_t pos = 0;
    size_t lastPos = 0;
    //std::cout << str << std::endl;
    while ( ( pos = str.find( delimiter, lastPos) ) != std::string::npos ){
        //std::cout << pos << " " << lastPos << " " << str.substr( lastPos, pos -lastPos) << std::endl;
        subStrings.push_back( str.substr( lastPos, pos - lastPos ) );
        lastPos = pos + 1;
    }
    //std::cout << pos << " " << lastPos << " " << str.substr( lastPos, pos -lastPos) << std::endl;
    subStrings.push_back( str.substr( lastPos) );
    return subStrings;
}

std::map < float, float > loadFloatMap( const std::string & filename ){
    std::map < float, float > aMap;

    std::fstream file; if ( !openInFile( filename, & file ) ){
        std::cerr << WHERE_AM_I << " Map not found: " << filename << std::endl;
        return aMap;
    }

    std::vector < std::string > row;
    while ( ! file.eof() ) {
        row = getNonEmptyRow( file );
        if ( row.size() == 2 ) aMap[ toFloat( row[ 0 ] ) ] = toFloat( row[ 1 ] );
    }

    file.close();
    return aMap;
}

std::map < int, int > loadIntMap( const std::string & filename ){
    std::map < int, int > aMap;

    std::fstream file; if ( !openInFile( filename, & file ) ){
        std::cerr << WHERE_AM_I << " Map not found: " << filename << std::endl;
        return aMap;
    }

    std::vector < std::string > row;
    while ( ! file.eof() ) {
        row = getNonEmptyRow( file );
        if ( row.size() == 2 ) aMap[ toInt( row[ 0 ] ) ] = toInt( row[ 1 ] );
    }

    file.close();
    return aMap;
}

} // namespace GIMLI

