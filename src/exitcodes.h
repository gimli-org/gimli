/***************************************************************************
 *   Copyright (C) 2007-2017 by the GIMLi development team       *
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
 
#ifndef _GIMLI_EXITCODES__H
#define _GIMLI_EXITCODES__H

// http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen/a00653.html
// here we will define some std::exception instead if exit codes

#include <stdexcept>

namespace GIMLI {
    
inline void throwToImplement(const std::string & errString){
#ifndef USE_EXIDCODES
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
#endif
}

inline void throwRangeError(int exitCode, const std::string & errString, int idx, int low, int high){
#ifndef USE_EXIDCODES
    std::stringstream str(errString);
    str << " " << idx << " [" << low << ".." << high << ")" << std::endl;
    throw std::out_of_range(str.str());
#else
    std::cerr << str.str() << std::endl;
    exit(exitCode);
#endif
}

inline void throwLengthError(int exitCode, const std::string & errString){
#ifndef USE_EXIDCODES
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
    exit(exitCode);
#endif
}

DLLEXPORT bool debug();

inline void throwError(int exitCode, const std::string & errString){
#ifndef USE_EXIDCODES
    if (!debug()){
        throw std::length_error(errString);
    } else {
        std::cerr << "Debug: " << errString << std::endl;
    }
    throw std::length_error(errString);
#else
    std::cerr << errString << std::endl;
    exit(exitCode);
#endif
}



// the following may be DEPRECATED

// Commandozeile ungültig, wird von getopt automatisch gesetzt
#define EXIT_DEFAULT 1

// Commandozeile zu kurz;
#define EXIT_COMANDLINE 2

// Datenfile kann nicht geöffnet werden
#define EXIT_DATACONTAINER_FILE 3

// Im Datenfile sind zu wenig Elektroden definiert
#define EXIT_DATACONTAINER_NELECS 10

// Im Datenfile ist die Spezifizierung der Elektroden ungültig (Elektroden Tokenliste)
#define EXIT_DATACONTAINER_ELECS_TOKEN 11

// Im Datenfile ist die Spezifizierung der Datentypen ungültig (Daten Tokenliste)
#define EXIT_DATACONTAINER_UNKNOWN_TOKEN 12

// Im Datenfile ist die Anzahl der Daten ungültig oder kann nicht gelesen werden
#define EXIT_DATACONTAINER_DATASIZE 13

// Im Datenfile gibt es ungültige Elektrodenkonfigurationen, Elektrodennummer ist ggf. größer als die Anzahl an Elektroden
#define EXIT_DATACONTAINER_ELECS 14

// Im Datenfile ist entweder u oder i nicht definiert oder gleich 0.0
#define EXIT_DATACONTAINER_NO_U_I 15

// Im Datenfile sind zu wenig Topographiepunkte definiert
#define EXIT_DATACONTAINER_NTOPO 16

// Im Datenfile ist keine Formatspezifikation
#define EXIT_DATACONTAINER_NO_DATAFORMAT 17

// Kann die Datei für den Export des Netzes und Models nicht erstellen
#define EXIT_MESH_EXPORT_FAILS 21

// Kann die InsertDatei für den Export des Netzes und Models nicht erstellen
#define EXIT_POLY_CIRC_INSERT_ADDITIONALS 22

// allgemeiner datei fehler, kann ein file entweder lesend oder schreibend nicht öffnen
#define EXIT_OPEN_FILE 22

// allgemeiner datei fehler, kann ein file entweder lesend oder schreibend nicht öffnen
#define EXIT_DCMODELL_NO_ANALYT 30

//** interne Fehler evntl Bugs, die entweder vorher entdeckt werden sollten und keine Bedienfehler sind


// noch nicht implementiert
#define EXIT_TO_IMPL 100

// Im Tinyvector ungueltiger indexzugriff
#define EXIT_TINYVECTOR_INDEX 101

// Im Tinyvector ungueltige Vektorgröße
#define EXIT_TINYVECTOR_SIZE 102

// Fehler im Netzcontainer
#define EXIT_MESH_NO_ELEMENT 103

// Fehler im FEMer, kein Attribut definiert
#define EXIT_FEM_NO_RHO 104

// dcmodelling can not find sources
#define EXIT_BERT_NO_SOURCES 105

// AMD Fehler
#define EXIT_AMD_INTERNAL 106

// Datacontainer size mismatch
#define EXIT_DATACONTAINER_SIZE 107

// Datamap size mismatch
#define EXIT_DATAMAP_SIZE 108

// Inversion no jacobian defined
#define EXIT_INVERSION_NO_J 109

// Mesh request for unknown node
#define EXIT_MESH_NO_NODE 110

// Mesh request for unknown boundary
#define EXIT_MESH_NO_BOUNDARY 111

// Mesh request for unknown cell
#define EXIT_MESH_NO_CELL 112

// sparsematrix pattern invalid
#define EXIT_SPARSE_INVALID 113

// sparsematrix size request invalid
#define EXIT_SPARSE_SIZE 114

// mesh 
#define EXIT_MESH_CELLTYPE_UNDEFINED 115

// mesh 
#define EXIT_MESH_BOUNDTYPE_UNDEFINED 116

// index error during sensmat calculation
#define EXIT_SENS_INDEX 117

// Mesh does not know any electrodes
#define EXIT_MESH_NO_ELECS 118

// Mesh does not know any electrodes
#define EXIT_MATRIX_SIZE_INVALID 119

// Mesh does not know any electrodes
#define EXIT_VECTOR_SIZE_INVALID 120

} // namespace GIMLI

#endif // _GIMLI_EXITCODES__H
