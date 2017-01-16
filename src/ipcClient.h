/***************************************************************************
 *   Copyright (C) 2006-2013 by the GIMLi development team       *
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

#ifndef _GIMLI_IPC_CLIENT__H
#define _GIMLI_IPC_CLIENT__H

#include "gimli.h"

#include <deque>

#ifndef USE_IPC
    #define USE_IPC 0
#elif USE_IPC
    #if HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP 
		#undef USE_IPC
        #define USE_IPC 0
		//#include <boost/interprocess/managed_shared_memory.hpp>
        //using namespace boost::interprocess;
    #else
        #define USE_IPC 0
    #endif
#endif

namespace GIMLI{

//! Inter process communication via shared memory.
/*! Inter process communication via shared memory.
 * For sample usage see: unittest/testGIMLiMisc.h:testIPCSHM()
 */
class IPCClientSHM{
public:
    enum TypeID { BOOL, INT, DOUBLE, ARRAY };

    IPCClientSHM( bool verbose = false )
        : verbose_( verbose ), initialized_( false ) {
    }

private:
    /*! Copyconstructor */
    IPCClientSHM( const IPCClientSHM & ipc ){
        THROW_TO_IMPL
    }

public:
    void setSegmentName( const std::string & segmentName ){
        name_ = segmentName;
#if USE_IPC
        segment_ = managed_shared_memory( open_or_create, segmentName.c_str(), 65536 );
        initialized_ = true;
#endif
    }

    template < class ValueType > void set( const std::string & name, ValueType val, TypeID type ){
        if ( !initialized_ ) return;
#if USE_IPC
        segment_.find_or_construct< TypeID >( (name + "-TypeID").c_str() )( type );

        switch ( type ){
            case BOOL:{
                bool * v = segment_.find_or_construct< bool >( name.c_str() )( );
                *v = (bool)val;
            } break;
            case INT:{
                int * v = segment_.find_or_construct< int >( name.c_str() )( );
                *v = (int)val;
            } break;
            case DOUBLE:{
                double * v = segment_.find_or_construct< double >( name.c_str() )( );
                *v = (double)val;
            } break;
            case ARRAY:{
//                  double * v = segment_.find_or_construct< double * >( name.c_str() )[ val.size() ]( );
// //                 *v = (double)val;
            } break;
        }
#endif
    }

    template < class ValueType > ValueType get( const std::string & name ){
       if ( !initialized_ ) return ValueType(0);
#if USE_IPC

        TypeID * t = segment_.find< TypeID >( (name + "-TypeID" ).c_str() ).first;
        if ( !t ){
            throwError(1, WHERE_AM_I + " No ipc value registered for name: " + name + "-TypeID" );
        }

        switch ( *t ){
            case BOOL:{
                bool * v = segment_.find< bool >( name.c_str() ).first;
                return ValueType( *v );
            } break;
            case INT:{
                int * v = segment_.find< int >( name.c_str() ).first;
                return ValueType( *v );
            } break;
            case DOUBLE:{
                double * v = segment_.find< double >( name.c_str() ).first;
                return ValueType( *v );
            } break;
            case ARRAY:{

            }break;
        }
#endif
        return ValueType(0);
    }

    inline void setBool( const std::string & name, bool val ){ set< bool >( name, val, BOOL ); }
    inline bool getBool( const std::string & name ){ return get< bool >( name ); }

    inline void setInt( const std::string &name, int val ){ set< int >( name, val, INT ); }
    inline int getInt( const std::string &name ){ return get< int>( name ); }

    inline void setDouble( const std::string &name, double val ){ set< double >( name, val, DOUBLE ); }
    inline double getDouble( const std::string &name ){ return get< double >( name ); }

    void info(){
        if ( !initialized_ ) {
            std::cout << "IPC(shared_memory) not initialized or supported" << std::endl;
            return;
        }
#if USE_IPC

        std::cout << "name: " << segment_.get_size() << std::endl;
        std::cout << "size: " << segment_.get_size() << std::endl;
        std::cout << "free: " << segment_.get_free_memory() << std::endl;
        std::cout << "sanity: " << segment_.check_sanity() << std::endl;
        std::cout << "nNamed: " << segment_.get_num_named_objects() << std::endl;

        typedef managed_shared_memory::const_named_iterator const_named_it;
        const_named_it named_beg = segment_.named_begin();
        const_named_it named_end = segment_.named_end();

        for(; named_beg != named_end; ++named_beg){
            //A pointer to the name of the named object
            const managed_shared_memory::char_type *name = named_beg->name();
            //The length of the name
            std::size_t name_len = named_beg->name_length();
            //A constant void pointer to the named object
            const void *value = named_beg->value();

            std::cout << "\t" << std::string( name, name_len)
                //<< " " << managed_shared_memory::get_instance_length( named_beg->value() )
                << " " << value << std::endl;
        }

        typedef managed_shared_memory::const_unique_iterator const_unique_it;
        const_unique_it unique_beg = segment_.unique_begin();
        const_unique_it unique_end = segment_.unique_end();

        std::cout << "nUnique" << segment_.get_num_unique_objects() << std::endl;
        for(; unique_beg != unique_end; ++unique_beg){
            //The typeid(T).name() of the unique object
            const char *typeid_name = unique_beg->name();
            //The length of the name
            std::size_t name_len = unique_beg->name_length();
            //A constant void pointer to the unique object
            const void *value = unique_beg->value();
            std::cout << "\t" << std::string( typeid_name, name_len) << " " << value << std::endl;
        }
#endif
    }

    void enlarge(){
        if ( !initialized_ ) return;
        THROW_TO_IMPL
//         managed_shared_memory::grow("MySharedMemory", 500);
//         managed_shared_memory shm(open_only, "MySharedMemory");
    }

    static void free( const std::string & name ){
#if USE_IPC
        std::cout << "freeing ipc workspace: " << name << std::endl;
        shared_memory_object::remove( name.c_str() );
#endif
    }

protected:
    bool verbose_;
    bool initialized_;

    std::string name_;
#if USE_IPC
    managed_shared_memory segment_;
#endif
};

// using boost::asio::ip::tcp;
//
// inline void writeToData( char * data, const std::string & val, int & count ){
//   memcpy( &data[ count ], & val.c_str()[0], sizeof( char ) * val.length() );
//   count += sizeof( char ) * val.length();
// }
// template < class T > void writeToData( char * data, const T & val, int & count ){
//   memcpy( &data[ count ], & val, sizeof( T ) );
//   count += sizeof( T );
// }
// template < class T > void readFromData( T & val, const char * data, int & count ){
//   memcpy( & val, & data[ count ], sizeof( T ) );
//   count += sizeof( T );
// }
//
//
// static const uint8 IPCValueTypeInt = 0;
// static const uint8 IPCValueTypeDouble = 1;
//
// /*! IPC Message
// data=uint16,uint8,char*,int8,ValueSize
// BodyLength,NameLength,Name*,IPCValType,Value
// */
// class IPCMessage{
// public:
//     enum { HeaderLength  = sizeof( uint16 ) };
//     enum { MaxBodyLength = 1024 };
//
//     IPCMessage( );
//
//     IPCMessage( const std::string & name, int val = 0 ){
//         init_( name, IPCValueTypeInt, (int64)val );
//     }
//
//     IPCMessage( const std::string & name, double val ){
//         init_( name, IPCValueTypeDouble, val );
//     }
//
//     /*! Construct message from valid char *. Used to parse received data. */
//     IPCMessage( const char * data );
//
//     /*! Copy contructor */
//     IPCMessage( const IPCMessage & msg );
//
//     /*! Default destructor */
//     ~IPCMessage( );
//
//     bool decodeHeader();
//
//     char * data() const { return data_; }
//     char * data() { return data_; }
//
//     char * body() const { return data_ + HeaderLength; }
//     char * body() { return data_ + HeaderLength; }
//
//     uint16 length() const { return HeaderLength + bodyLength_; }
//     uint16 bodyLength() const { return bodyLength_; }
//
// protected:
//
//     template < class ValueType > void init_( const std::string & name, uint8 valType, const ValueType & val ){
//         size_t valSize = 1;
//         switch ( valType ){
//             case IPCValueTypeInt:   valSize = sizeof( int64 ); break;
//             case IPCValueTypeDouble: valSize = sizeof( double ); break;
//         }
//
//         //std::cout << "Out: " << name <<  " " << int( valType )<< " " << val << std::endl;
//         bodyLength_ = sizeof( uint8 ) + name.length() + sizeof( uint8 ) + valSize;
//
//         alloc_( bodyLength_ + HeaderLength );
//         int count = 0;
//
//         writeToData( data_, bodyLength_, count );
//         writeToData( data_, (uint8)name.length(), count );
//         writeToData( data_, name, count );
//         writeToData( data_, valType, count );
//         writeToData( data_, val, count );
//
// //         std::cout << "#";
// //         for ( size_t i = 0; i < length(); i ++ ) std::cout << (unsigned int) data_[ i ] << " ";
// //         std::cout << "#" << std::endl;
//     }
//
//     void alloc_( size_t l ){
//         data_ = new char[ HeaderLength + l ];
//     }
//
//     char * data_;
//     uint16 bodyLength_;
//
//     uint8 valueType_;
//     std::string name_;
// };
//
// /*! Client for inter process communication (IPC). Send messages to a server using sockets.
//     Very basic messages allowed at the moment, you can just send or receive namend values to/from the server
// */
// class IPCClient{
// public:
//
//     typedef std::deque< IPCMessage > MessageQueue;
//
//     /*! Create an IPCClient */
//     IPCClient( bool verbose = false );
//
//     ~IPCClient();
//
//     /*! Connect to the server host:port.*/
//     void connect( const std::string & hostName, int port );
//
//     /*! Send value val of type (int or double) with a given name (maximal name length =256chars) to the server. */
//     template < class ValueType > void send( const std::string & name, const ValueType & val ){
//         this->handleSendMsg_( IPCMessage( name, val ) );
//     }
//
//     void setVerbose( bool verbose ) { verbose_ = verbose; }
//
//     bool verbose( ) const { return verbose_; }
//
// protected:
//     void handleConnect_( const boost::system::error_code & e, boost::asio::ip::tcp::resolver::iterator endpoint_iterator );
//     void handleSendMsg_( const IPCMessage & body );
//     void handleWrite_( const boost::system::error_code & e );
//
//     void handleReadBody(const boost::system::error_code & e );
//     void handleReadHeader(const boost::system::error_code & e );
//
//     void closeConnection_();
//
//     boost::asio::io_service & io_service_;
//     tcp::socket socket_;
//
//     bool verbose_;
//     bool online_;
//
//     MessageQueue outMessageQueue_;
//     IPCMessage incomingMsg_;
// };



} // namespace GIMLI{

#endif //_GIMLI_IPC_CLIENT__H
