/***************************************************************************
 *   Copyright (C) 2006-2008 by the resistivity.net development team       *
 *   Carsten Ruecker carsten@resistivity.net                                *
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

#include "platform.h"

#include <wx/window.h>

#ifdef __WXGTK__
    // Needed for random GTK/GDK macros and functions
    #include <gdk/gdkx.h>
    #include <gtk/gtk.h>
#endif

#include <sstream> 
namespace GIMLI{
std::string getGTKWindowHandleStr( void * window ) {
    std::string handle;
    if ( window ) {
            // format display:screen:window, which has variable types ulong:uint:ulong.
            GtkWidget* widget = reinterpret_cast<wxWindow*>(window)->GetHandle();
            gtk_widget_realize( widget );   // Mandatory. Otherwise, a segfault happens.
            //Display* display = GDK_WINDOW_XDISPLAY( widget->window );
            Window wid = GDK_WINDOW_XWINDOW( widget->window );      // Window is a typedef for XID, which is a typedef for unsigned int
            /* Get the right display (DisplayString() returns ":display.screen") */
            //std::string displayStr = DisplayString( display );
            //displayStr = displayStr.substr( 1, ( displayStr.find( ".", 0 ) - 1 ) );
            /* Put all together */
            std::stringstream handleStream;
            handleStream << 0 << ':' << 0 << ':' << wid;
            //handleStream << displayStr << ':' << DefaultScreen( display ) << ':' << wid;
            handle = handleStream.str();
            return handle;
    } else {
        return "0";
    }
}
}