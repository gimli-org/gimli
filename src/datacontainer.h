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

template < class T > bool checkValue( T val, T min , T max, T tol ){
    return !( std::fabs( val ) < tol || val < min || val > max || std::isinf( val ) || std::isnan( val ) );
}

class DLLEXPORT DataContainer{
public:
    /*! Simple Constructor, builds an empty data container.*/
    DataContainer( );

    /*! Constructor, builds a data container, and fills the data from the file fileName. See \ref load.*/
    DataContainer( const std::string & fileName );

    /*! Copy Constructor */
    DataContainer( const DataContainer & data );

    /*! Destructor. */
    virtual ~DataContainer( );

    /*! Assigns a copy of data and returns a reference to this DataContainer */
    DataContainer & operator = ( const DataContainer & data );

    /*! Init the datacontainer for dc-electrical needs. (inti a, b ,m ,b and valid ) */
    virtual void init();

    /*! Init a token translator map. store a map ['alias'] = 'token'.
        Only lowercase starting aliases are handled i.e. alias( power, u ) map to token => u
        On request the first char of the alias will converted to lowercase. e.g.
        translateToken( 'power' ) -> return 'u'
        translateToken( 'Power' ) -> return 'u'
        translateToken( 'u' ) -> return 'u'
        translateToken( 'U' ) -> return 'u' */
    virtual void initTokenTranslator();

    /*! Return true if there is a valid translation for the token */
    bool haveTranslationForAlias( const std::string & alias ) const;

    /*! Return the token translation */
    std::string translateAlias( const std::string & alias ) const;

    /*! Clear the container, remove all electrodes and data. */
    virtual void clear();

    inline const RVector & operator() ( const std::string & token ) const { return get( token ); }

    inline RVector & operator() ( const std::string & token ) { return *ref( token ); }

    inline size_t size( ) const { return dataMap_.find( "valid" )->second.size(); }

    void createFourPointData( uint i, int a, int b, int m, int n );

    inline const std::vector < RVector3 > & additionalPoints() const { return topoPoints_; }

    inline const std::map< std::string, RVector > & dataMap() const { return dataMap_; }

    /*! Return true if exists and all elements != 0.0 */
    inline bool nonZero( const std::string & token ) const {
        if ( exists( token ) ) return ( min( abs( dataMap_.find( token )->second ) ) > TOLERANCE );
        return false;
    }

    /*! Return true if exists and at least one value != 0.0 */
    inline bool haveData( const std::string & token ) const {
        if ( exists( token ) ) return ( max( abs( dataMap_.find( token )->second ) ) > TOLERANCE );
        return false;
    }

    /*! Return true if the data with the token exists. */
    inline bool exists( const std::string & token ) const { return dataMap_.count( token ); }

    inline const std::string & inputFormatString() const { return inputFormatString_; }

    /*! Return reference to the token translator map */
    inline std::map< std::string, std::string > & tokenTranslator() { return tT_; }

    /*! Loads the data from a file. See. http://www.resistivity.net/?unidata for details on the fileformat.*/
    virtual int load( const std::string & fileName );

    /*! Save the data to the file fileName. formatVals -- data tokenlist to save; formatElecs -- electrodes tokenlist to save 
     Saves only valid data. */
    virtual int save( const std::string & fileName, const std::string & formatVals, const std::string & formatElecs, bool verbose = false ) const;

    /*! Shortcut for save.*/
    inline int save( const std::string & fileName, const std::string & formatVals, bool verbose = false ) const {
        return save( fileName, formatVals, "x y z", verbose ); }

    /*! Shortcut for save */
    inline int save( const std::string & fileName, bool verbose = false ) const {
        return save( fileName, "a b m n rhoa err k", verbose ); }

    /*! Shows some information belongs to the DataContainer.*/
    void showInfos() const ;

    void resize( uint size );

    /*! Return string with space separated list of all available tokens. 
     */
    std::string tokenList() const;

    /*! Add new data field with optional description. Throws an exception if the data field size is not the same size of the DataContainer.
        \param token String to identify the data
        \param data \ref RVector of the data
        \param description String that describe the data.
    */
    void add( const std::string & token, const RVector & data, const std::string & description = "" );

    /*! Set the data for a given token. If there is no such data, new data will be added. 
     * Throws an exception if the data field size is not the same size of the \ref DataContainer.
        \param token String to identify the data
        \param data \ref RVector of the data
    */
    void set( const std::string & token, const RVector & data );

    /*! Read only access to a data field. Throws an exception if data doesn't exist.
     * \param token String to identify the data
    */
    const RVector & get( const std::string & token ) const ;

    /*! Read/Write access via a pointer to the data field
     * \param token String to identify the data
    */
    RVector * ref( const std::string & token );

    /*! Set description string for a specific data field. If the data doesn't exist nothing is done.
        \param token String that identify the data to be described.
        \param description String that describe the data.
    */
    void setDataDescription( const std::string & token, const std::string & description );
    
    /*! Returns a copy of the description string for the specified data field. Return empty string if the data doesn't exist.
        \param token String that identify the data.
    */
    std::string dataDescription( const std::string & token ) const;
    
    /*! Inline remove data from index vector. Remove all data that are covered by idx. Electrodes are preserved. */
    void remove( const std::vector < size_t > & idx );


    /*! Mark data valid by index vector. Shortcut for this->ref("valid")->setVal( idx, valid ). */
    inline void markValid( const std::vector < size_t > & idx, bool valid = true ){
        dataMap_[ "valid" ].setVal( valid, idx );
    }

    /*! Mark single data valid. this->ref("valid")->setVal( idx, valid ). */
    inline void markValid( size_t idx, bool valid = true ){
        dataMap_[ "valid" ].setVal( valid, idx );
    }

    template < class T > void remove( const std::string & token, T min, T max, T tol, bool markOnly = false){

        //this->markValid( find( this->get( token ) < min && )
        //** still little hackish will be better with dataMap_;

        for ( uint i = 0; i < this->size(); i ++ ) {
            if ( !checkValue( get(token)[ i ], min, max, tol ) ) {
                dataMap_["valid"][ i ] = false;
            }
        }
        if ( ! markOnly ) removeInvalid();
    }

    void removeInvalid();

    /*! Create new DataContainer from index vector. Create a copy and remove all data that are not covered by idx. Electrodes are preserved. */
    DataContainer filter( const std::vector < size_t > & idx ) const ;

protected:
    virtual void copy_( const DataContainer & data );

    std::vector < RVector3 >         topoPoints_;

    std::string inputFormatString_;

    //! Data fields
    /*! Store the data field map < token, data >. All data fields have the same size. */
    std::map< std::string, RVector > dataMap_; 
    
    //! Description for the data fields
    /*! Store an optional description for the associated data field < token, description >. \warning Description will not yet be saved. */
    std::map< std::string, std::string > dataDescription_;  
    
    //** tokenTranslater for renaming formats to known cases:
    std::map< std::string, std::string > tT_;

}; // class DataContainer

} // namespace GIMLI{

#endif // DATACONTAINER__H
