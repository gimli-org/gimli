/***************************************************************************
 *   Copyright (C) 2006-2011 by the resistivity.net development team       *
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

#ifndef _GIMLI_DATACONTAINER__H
#define _GIMLI_DATACONTAINER__H

#include "gimli.h"
#include "baseentity.h"
#include "vector.h"

#include <vector>
#include <set>
#include <map>

namespace GIMLI{

//! DataContainer to store, load and save data in the gimli unified data format
/*! DataContainer to store, load and save data in the gimli unified data format. 
 The DataContainer contains a data map that holds the data itself. Each map entry can be identified by tokens.
 By default there is a data field with the token 'valid' to mark the validity of the data.
 There is also a vector of unique sensor positions that holds sensor locations and a set of additional points, e.g., topographic positions. 
 A vector of indices to the sensor positions can be defined for each datum. e.g., Store an index-vector 'a' to the first current electrode 'A' of a ERT measurement. 
 If you need a special DataContainer you should derive a child and specify token translator and sensor index entries.
 There is also a unit test for the DataContainer that may help to understand what is it good for. */
class DLLEXPORT DataContainer{
public:
    /*! Simple Constructor, builds an empty data container.*/
    DataContainer( );

    /*! Constructor, builds a data container and fills the data from a file. See \ref load.
      \param fileName String of the file name */
    DataContainer( const std::string & fileName );

    /*! Copy Constructor */
    DataContainer( const DataContainer & data );

    /*! Destructor. */
    virtual ~DataContainer( );

    /*! Assign a copy of data and return a reference to this DataContainer */
    DataContainer & operator = ( const DataContainer & data );

    /*! Return read only reference to the RVector at the data map that is associated to the token */
    inline const RVector & operator() ( const std::string & token ) const { return get( token ); }

    /*! Return reference to the RVector at the data map that is assiciated to the token */
    inline RVector & operator() ( const std::string & token ) { return *ref( token ); }
    
    /*! Init the datacontainer for dc-electrical needs. (int a, b ,m ,b and valid ) */
    virtual void init();

    /*! Init a token translator map. store a map ['alias'] = 'token'.
        Only lowercase starting aliases are handled i.e. alias( power, u ) map to token => u
        On request, the first char of the alias will converted to lowercase. e.g.
        translateToken( 'power' ) -> return 'u'
        translateToken( 'Power' ) -> return 'u'
        translateToken( 'u' ) -> return 'u'
        translateToken( 'U' ) -> return 'u' */
    virtual void initTokenTranslator();

    /*! Return true if there is a valid translation for the token */
    bool haveTranslationForAlias( const std::string & alias ) const;

    /*! Return the token translation */
    std::string translateAlias( const std::string & alias ) const;

    /*! Clear the container, remove all sensor locations and data. */
    virtual void clear();

    /*! Return the size of the data map */
    inline size_t size( ) const { return dataMap_.find( "valid" )->second.size(); }

    /*! Return the complete data map as read-only map */
    inline const std::map< std::string, RVector > & dataMap() const { return dataMap_; }
    
    // START Sensor related section
    /*! Return the complete sensor positions as read-only */
    inline const std::vector< RVector3 > & sensorPositions() const { return sensorPoints_; }

    /*! Return a single sensor position */
    inline const RVector3 & sensorPosition( uint i ) const { return sensorPoints_[i]; }

    /*! Create a valid sensor at a given position and returns the id of the sensor. 
        Is there already a sensor at the given position NO new sensor will be created. 
        Atm. brute force search with a snapping distance of tolerance is done.
        \param pos RVector3 of the sensor position
        \param tolerance Double of the snapping tolerance */
    long createSensor( const RVector3 & pos, double tolerance = 1e-3 );
    
    /*! Return the number of sensors. */
    uint sensorCount(  ) const { return sensorPoints_.size(); }
 
    /*! Mark the data field entry as sensor index. */
    void registerSensorIndex( const std::string & token );
    
    /*! Return true if the field entry is of type sensor index. */
    bool isSensorIndex( const std::string & token ) const ;

    /*! Return true if the sensor indices on a loaded/saved file starting 1. Internally the indices stored from 0. */
    bool sensorIndexOnFileFromOne() const { return sensorIndexOnFileFromOne_ ;}
        
    /*! Mark all data invalid that use a sensor index greater than sensor count. */
    void markInvalidSensorIndices(); 
    
    /*! Remove all data that contains the sensor and the sensor itself. 
    *\param idx uint idx single index for a sensor regarding sensorPoints_
    */
    void removeSensorIdx( uint idx );
    
    /*! Remove all data that contains the sensor and the sensor itself. *
     *\param idx IndexArray array of indices regarding sensorPoints_
     */
    void removeSensorIdx( const IndexArray & idx );
    
    // END Sensor related section
    
    /*! Return the additional points */
    inline const std::vector < RVector3 > & additionalPoints() const { return topoPoints_; }
    
    /*! Return true if token data exist and all elements != 0.0 */
    inline bool nonZero( const std::string & token ) const {
        if ( exists( token ) ) return ( min( abs( dataMap_.find( token )->second ) ) > TOLERANCE );
        return false;
    }

    /*! Return true if token data exist and at least one value != 0.0 */
    inline bool haveData( const std::string & token ) const {
        if ( exists( token ) ) return ( max( abs( dataMap_.find( token )->second ) ) > TOLERANCE );
        return false;
    }

    /*! Return true if the data with the token exists. */
    inline bool exists( const std::string & token ) const { return dataMap_.count( token ); }

    /*! Return the token list of a previously loaded file. */
    inline const std::string & inputFormatString() const { return inputFormatString_; }

    /*! Return reference to the token translator map. */
    inline std::map< std::string, std::string > & tokenTranslator() { return tT_; }

    /*! Loads the data from a file. See save for details on the fileformat.*/
    virtual int load( const std::string & fileName );

    /*! Save the data to a file. Saves only valid data(except formatData == "all"). File format is\n\n
     * Number of Sensors\n
     * #Sensor tokens\n
     * Sensor[token 1] Sensor[token n]\n
     * Number of Data\n
     * #Data tokens\n
     * Data[token 1] Data [token n]\n
     * Number of additional points\n
     * Additional points\n
     * 
     * http://www.resistivity.net/?unidata\n
        \param fileName String of the file name 
        \param formatData String to specify the tokens of the data map to be save. If formatData == "all", all datafields will be saved inclusive invalid data.
        \param formatSensor String to specify the tokens of the sensor format to be save */
    virtual int save( const std::string & fileName, const std::string & formatData, const std::string & formatSensor, bool verbose = false ) const;

    /*! Shortcut for \ref save( fileName, formatData, "x y z", verbose ); */
    inline int save( const std::string & fileName, const std::string & formatData = "all", bool verbose = false ) const {
        return save( fileName, formatData, "x y z", verbose ); }

    /*! Show some information that belongs to the DataContainer.*/
    void showInfos() const ;

    /*! Resize the data map and all data fields to size.*/
    void resize( uint size );

    /*! Return string with space separated list of all available data tokens. */
    std::string tokenList() const;

    /*! Add new data field with optional description. Throws an exception if the data field size is not the same size of the DataContainer.
        \param token String to identify the data
        \param data \ref RVector of the data
        \param description String that describe the data */
    void add( const std::string & token, const RVector & data, const std::string & description = "" );

    /*! Set the data for a given token. If there is no such data, new data will be added. 
     * Throws an exception if the data field size is not the same size of the \ref DataContainer.
        \param token String to identify the data
        \param data \ref RVector of the data */
    void set( const std::string & token, const RVector & data );

    /*! Read only access to a data field. Throws an exception if token data don't exist.
     * \param token String to identify the data */
    const RVector & get( const std::string & token ) const ;

    /*! Read/Write access via a pointer to the token data field
     * \param token String to identify the data */
    RVector * ref( const std::string & token );

    /*! Set description string for a specific data field. If the data doesn't exist nothing is done.
        \param token String that identify the data to be described
        \param description String that describe the data */
    void setDataDescription( const std::string & token, const std::string & description );
    
    /*! Returns a copy of the description string for the specified data field. Return empty string if the data doesn't exist.
        \param token String that identify the data.
    */
    std::string dataDescription( const std::string & token ) const;
    
    /*! Inplace remove data from index vector. Remove all data that are covered by idx. Sensors are preserved.*/
    virtual void remove( const std::vector < size_t > & idx );

    /*! Create new DataContainer that only contains the values that are covered by idx. Sensors are preserved.*/
    virtual DataContainer filter( const std::vector < size_t > & idx ) const ;
    
    /*! Mark data valid by index vector. Shortcut for this->ref("valid")->setVal( idx, valid ). */
    inline void markValid( const std::vector < size_t > & idx, bool valid = true ){
        dataMap_[ "valid" ].setVal( valid, idx );
    }
    
    /*! Mark single data valid. this->ref("valid")->setVal( idx, valid ). */
    inline void markValid( size_t idx, bool valid = true ){
        dataMap_[ "valid" ].setVal( valid, idx );
    }

    /*! Remove all data[valid] == 0. Sensors are preserved.*/
    void removeInvalid();

    /*! Remove all unused sensors from this DataContainer and recount data sensor index entries */
    void removeUnusedSensors();

protected:
    virtual void copy_( const DataContainer & data );

    std::string inputFormatString_;

    //! Data map
    /*! Store the data map < token, data >. All data map entries have the same size. */
    std::map< std::string, RVector > dataMap_; 
    
    //! Sensor positions
    /*! Store the sensor positions. */
    std::vector < RVector3 > sensorPoints_; 

    //! Data field that is sensor index
    /*! Store the field token that represent sensor indices */
    std::set< std::string > dataSensorIdx_; 
    
    //! Description for the data map entries
    /*! Store an optional description for the associated data field < token, description >. \warning Description will not yet be saved. */
    std::map< std::string, std::string > dataDescription_;  
    
    /*! Store additionally points */
    std::vector < RVector3 > topoPoints_;

    /*! tokenTranslator for renaming formats to known cases */
    std::map< std::string, std::string > tT_;
    
    /*! Determine if the sensor indices should start from 0 or 1 */
    bool sensorIndexOnFileFromOne_;

}; // class DataContainer

} // namespace GIMLI

#endif // DATACONTAINER__H
