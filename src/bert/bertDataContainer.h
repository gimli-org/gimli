/***************************************************************************
 *   Copyright (C) 2011-2017 by the resistivity.net development team       *
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

#ifndef _BERT_BERT_DATACONTAINER__H
#define _BERT_BERT_DATACONTAINER__H

#include "bert.h"

#include <datacontainer.h>

namespace GIMLI{

typedef std::pair< int, int > CurrentPattern;

class DLLEXPORT DataContainerERT : public GIMLI::DataContainer{
public:
    /*! Simple Constructor, builds an empty data container.*/
    DataContainerERT();

    /*! Constructor, builds a data container and fills the data from a file. See \ref load.
      \param fileName String of the file name */
    DataContainerERT(const std::string & fileName, bool removeInvalid=true);

    /*! Copy Constructor */
    DataContainerERT(const DataContainerERT & data);

    /*!
     * Default destructor
     */
    virtual ~DataContainerERT();

    /*!
     * Define ert sensorindex names 'a b m n' and ert data fields 'u i r rhoa err k ip iperr'
     */
    virtual void init();

    /*!
     * Define ert token translators
     */
    virtual void initTokenTranslator();

    void createFourPointData(Index i, long eaID, long ebID, long emID, long enID);

    virtual void checkDataValidityLocal();

    CurrentPattern currentPatternToElectrode(SIndex pattern);

    SIndex electrodeToCurrentPattern(SIndex a, SIndex b) const;

    std::set < SIndex > currentPattern(bool reciprocity=false);

    /*! Merge duplicate data by averaging. Sort the DataContainerERT as well.*/
    void averageDuplicateData(bool verbose=false);

    void fitFillSize();

protected:

    Index fillCounter_;
};

} // namespace BERT

#endif // _BERT_BERT_DATACONTAINER__H

