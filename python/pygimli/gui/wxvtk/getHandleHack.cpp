// Copyright (C) 2007 Maryland Robotics Club 
// Copyright (C) 2007 Joseph Lisee <jlisee@umd.edu> 
// All rights reserved. 
// 
// Author: Joseph Lisee <jlisee@umd.edu> 
// File:  utility/wxogre_util/src/wxogre_util.h 

// modified for pybert's vtk import to fit python-2.9@linux

#include <sstream> 
#include <string> 
#include <cstring> 
#include <iostream> 
using namespace std; 

// wxWidgets includes 
#include <wx/window.h> 

// Platform specific includes 
#ifdef __WXGTK__ 
// Needed for random GTK/GDK macros and functions 
#include <gdk/gdkx.h> 
#include <gtk/gtk.h> 

#endif


extern "C" 
{ 
  /** Returns a window handle string from a wxWindow in Ogre form 
     @param window 
         Must be a pointer to a wxWindow or its subclass. 
     @return 
        The window handle in its display:xid form. 
   */ 
  extern long get_window_handle(void* window);
} 

long get_window_handle( void * window ) { 
    std::string handle;
    if ( window ) { 
        GtkWidget* widget = reinterpret_cast<wxWindow*>(window)->m_widget;//->GetHandle();
        if ( !widget ){
            return -1;
        }

        gtk_widget_realize( widget );   // Mandatory. Otherwise, a segfault happens.
        Window wid = GDK_WINDOW_XWINDOW( widget->window );      // Window is a typedef for XID, which is a typedef for unsigned int

        return (long)wid;
    } else {
        return -2;
    }
    
    return 0;
} 

