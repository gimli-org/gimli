/******************************************************************************
 *   Copyright (C) 2006-2019 by the GIMLi development team                    *
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

#include "ipcClient.h"

//#include <boost/thread.hpp>
#include <iostream>

namespace GIMLI{

// boost::asio::io_service io_service;
// boost::mutex __sendMsgMutex__;
//
// IPCMessage::IPCMessage( ){
//     bodyLength_ = MaxBodyLength;
//     alloc_( bodyLength_ + HeaderLength );
// }
//
// IPCMessage::IPCMessage( const char * data ) : data_( NULL ){
//     int count = 0;
//
//     uint8 nameLength; readFromData( nameLength, data, count );
//     std::string msg( (char*) &data[ count ], nameLength ); count += nameLength;
//     name_ = msg;
//     readFromData( valueType_, data, count );
//     std::cout << "In: " << name_ << " " << int( valueType_ ) << " ";
//
//     switch ( valueType_ ){
//         case IPCValueTypeInt:{
//             int64 val; readFromData( val, data, count );
//             std::cout << val << std::endl;
//         } break;
//         case IPCValueTypeDouble: {
//             double val; readFromData( val, data, count );
//             std::cout << val << std::endl;
//         } break;
//         break;
//     }
// }
//
// IPCMessage::IPCMessage( const IPCMessage & msg ){
//     bodyLength_ = msg.bodyLength();
//     alloc_( msg.length() );
//     memcpy( data_, msg.data(), msg.length() );
// }
//
// IPCMessage::~IPCMessage( ){
//     if ( data_ ) delete [] data_;
// }
//
// bool IPCMessage::decodeHeader(){
//     int count = 0;
//     readFromData( bodyLength_, data_, count );
//     if ( bodyLength_ > MaxBodyLength || bodyLength_ == 0) {
//         std::cerr << "Message::decodeHeader fails: " << bodyLength_ << " > " << MaxBodyLength << std::endl;
//         bodyLength_ = 0;
//         return false;
//     }
//     return true;
// }
//
// IPCClient::IPCClient( bool verbose )
// : io_service_( io_service ), socket_( io_service ), verbose_( verbose ), online_( false ) {
// }
//
// IPCClient::~IPCClient(  ){
//     this->closeConnection_();
// }
//
// void IPCClient::connect( const std::string & hostName, int port ){
//     if ( verbose_ ) std::cout << "Connecting to " << hostName  << ":" << toStr( (int)port ) << std::endl;
//     try{
//         boost::asio::ip::tcp::resolver resolver( io_service_ );
//         boost::asio::ip::tcp::resolver::query query( hostName, toStr( (int)port) );
//
//         boost::asio::ip::tcp::resolver::iterator endpoint_iterator = resolver.resolve( query );
//         boost::asio::ip::tcp::endpoint endpoint = *endpoint_iterator;
//         online_ = true;
//         socket_.async_connect( endpoint, boost::bind( &IPCClient::handleConnect_, this, boost::asio::placeholders::error, ++endpoint_iterator ) );
//     } catch ( std::exception & e) {
//         std::cerr << std::string( "std::exception: " ) << e.what() << std::endl;
//     }
//
//     if ( online_ ){
//         boost::thread t( boost::bind( & boost::asio::io_service::run, & io_service_ ) );
//     }
// }
//
// void IPCClient::handleConnect_( const boost::system::error_code & e, boost::asio::ip::tcp::resolver::iterator endpoint_iterator ){
//     if ( !e ) {
//         if ( verbose_ ) std::cout << "Connected: waiting for welcome. " << std::endl;
//         //std::cout << "IPCMessage::readHeader" << std::endl;
//         boost::asio::async_read( socket_,
//                                 boost::asio::buffer( incomingMsg_.data(), IPCMessage::HeaderLength ),
//                                 boost::bind( & IPCClient::handleReadHeader, this, boost::asio::placeholders::error ) );
//     }
//     else if ( endpoint_iterator != boost::asio::ip::tcp::resolver::iterator() ) {
//         // Try the next endpoint.
//         socket_.close();
//         boost::asio::ip::tcp::endpoint endpoint = *endpoint_iterator;
//         socket_.async_connect( endpoint, boost::bind( &IPCClient::handleConnect_, this, boost::asio::placeholders::error, ++endpoint_iterator ) );
//     } else {
//         std::cerr << "IPCClient::handleConnect: " << e.message() << std::endl;
//     }
// }
//
// void IPCClient::handleReadHeader(const boost::system::error_code & e){
//     if ( !e && incomingMsg_.decodeHeader() ){
//         //std::cout << "IPCMessage::readBody " << incomingMsg_.bodyLength() << std::endl;
//         boost::asio::async_read( socket_,
//                                     boost::asio::buffer( incomingMsg_.body(), incomingMsg_.bodyLength() ),
//                                     boost::bind( & IPCClient::handleReadBody, this, boost::asio::placeholders::error));
//     } else {
//         std::cerr << "IPCClient::handleReadHeader() " << e.message() << std::endl;
//         this->closeConnection_();
//     }
// }
//
// void IPCClient::handleReadBody(const boost::system::error_code & e ){
//     if ( !e ){
// //        std::cout << incomingMsg_.bodyLength() << std::endl;
//         //std::cout << "IPCClient::handleReadBody()" << std::endl;
//         IPCMessage( incomingMsg_.body() );
//         //std::cout << "IPCMessage::readHeader" << std::endl;
//         boost::asio::async_read( socket_,
//                                     boost::asio::buffer( incomingMsg_.data(), incomingMsg_.HeaderLength ),
//                                     boost::bind( & IPCClient::handleReadHeader, this, boost::asio::placeholders::error ) );
//     } else {
//         std::cerr << "IPCClient::handleReadBody() " << e.message() << std::endl;
//         this->closeConnection_();
//     }
// }
//
// void IPCClient::handleSendMsg_( const IPCMessage & msg ) {
//     boost::mutex::scoped_lock lock( __sendMsgMutex__ );
//     if ( online_ ){
//         bool write_in_progress = !outMessageQueue_.empty();
//         //std::cout << "push_back: " << (size_t)msg.body()[0] << " " << write_in_progress << " " << outMessageQueue_.size() << std::endl;
//         {
//             //boost::mutex::scoped_lock lock( __sendMsgMutex__ );
//             outMessageQueue_.push_back( IPCMessage( msg ) );
//         }
//
//         if ( !write_in_progress ) {
//             //std::cout << "write1: " << outMessageQueue_.front().length() << std::endl;
// //             boost::asio::async_write( socket_,
// //                                         boost::asio::buffer( outMessageQueue_.front().data(), outMessageQueue_.front().length() ),
// //                                         boost::bind( &IPCClient::handleWrite_, this, boost::asio::placeholders::error ) );
//             boost::asio::write( socket_,boost::asio::buffer( outMessageQueue_.front().data(), outMessageQueue_.front().length() ) );
//             handleWrite_( boost::system::error_code() );
//
//         }
//     }
// }
//
// void IPCClient::handleWrite_( const boost::system::error_code & e ) {
//     //std::cout << "handleWrite_1: " << outMessageQueue_.size() << std::endl;
//     if ( online_ ){
//       //  std::cout << "handleWrite_2: " << outMessageQueue_.size() << std::endl;
//         if ( !e ){
//             //! write until queue is empty
//             if ( !outMessageQueue_.empty() ) {
//               //  boost::mutex::scoped_lock lock( __sendMsgMutex__ );
//                 //std::cout << "pop_front: " << (size_t)outMessageQueue_.front().body()[0] << std::endl;
//                 outMessageQueue_.pop_front();
//             }
//             if ( !outMessageQueue_.empty() ) {
//                 //std::cout << "write2: " << (size_t)outMessageQueue_.front().body()[0]<< std::endl;
// //                  boost::asio::async_write( socket_,
// //                                             boost::asio::buffer( outMessageQueue_.front().data(), outMessageQueue_.front().length() ),
// //                                             boost::bind( &IPCClient::handleWrite_, this, boost::asio::placeholders::error ) );
//                 boost::asio::write( socket_,boost::asio::buffer( outMessageQueue_.front().data(), outMessageQueue_.front().length() ) );
//                 handleWrite_( boost::system::error_code() );
//             }
//
//         } else {
//             std::cerr << "IPCClient::handleWrite: " << e.message() << std::endl;
//             closeConnection_();
//         }
//     }
// }
//
// void IPCClient::closeConnection_() {
//     if ( online_ ){
//         try {
//             socket_.close();
//         } catch ( std::exception & e ){
//             std::cerr << "IPCClient::closeConnection_() " << e.what() << std::endl;
//         }
//         online_ = false;
//     }
// }

} // namespace GIMLI
