
/******************************************************************************
 *   Copyright (C) 2006-2024 by the GIMLi development team                    *
 *   Carsten Rücker carsten@gimli.org                                         *
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
#pragma once

#include "gimli.h"

namespace GIMLI{

/*! Simple hexdump.
borrowed from from https://stackoverflow.com/questions/29242/off-the-shelf-c-hex-dump-code
*/
DLLEXPORT void hexdump(void *ptr, Index buflen);

class ByteBuffer{
public:
    ByteBuffer(){

    }
    void push_back(unsigned char c){
        _buf.push_back(c);
    }
    inline char * begin() const { return (char *)&_buf[0]; }
    inline char * end() const { return begin() + size(); }
    inline Index size() const { return _buf.size(); }

#ifndef PYGIMLI_CAST
    std::vector < uint8 > & buf() { return _buf; }
#endif
protected:
    std::vector < uint8 > _buf;

};

}