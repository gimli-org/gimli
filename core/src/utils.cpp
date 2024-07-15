
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

#include "utils.h"


namespace GIMLI{

void hexdump(void *ptr, Index buflen) {
    unsigned char *buf = (unsigned char*)ptr;

    for (Index i = 0; i < buflen; i+=16) {
        printf("%06x: ", (int)i);
        for (Index j = 0; j < 16; j++){
            if (i + j < buflen)
                printf("%02x ", buf[i+j]);
            else
                printf("   ");
        }
        printf(" ");
        for (Index j = 0; j < 16; j++){
            if (i + j < buflen)
                printf("%c", isprint(buf[i+j]) ? buf[i+j] : '.');
        }
        printf("\n");
    }
}

}