# -*- coding: utf-8 -*-

import sys
import os

try:
    import wx
    import wx.xrc as xrc
except ImportError as e:
    import traceback
    #traceback.print_exc(file=sys.stdout)
    sys.stderr.write("No proper wx installed'.\n")
try:
    from agw import aui
    from agw.aui import aui_switcherdialog as ASD
except ImportError: # if it's not there locally, try the wxPython lib.

    try:
        import wx.lib.agw.aui as aui
        from wx.lib.agw.aui import aui_switcherdialog as ASD
    except ImportError: # if it's not there locally, try the wxPython lib.
        sys.stderr.write("No proper wx.aui installed'.\n")

from pygimli.utils import unicodeToAscii

class Property:
    """What is this??"""
    def __init__( self, name, default = None, valType = None, ctrl = None, ctrlEvent = None, targetFunct = None ):
        self.name        = name
        self.val         = default
        self.valType     = valType

        self.ctrl        = ctrl
        self.ctrlEvent   = ctrlEvent
        self.targetFunct = targetFunct

        self.minVal      = None
        self.maxVal      = None

        self.parent      = None

    def __call__( self ):
        return self.val

    def setValidRange( self, minVal, maxVal ):
        self.minVal = minVal
        self.maxVal = maxVal

    def setCtrl( self, ctrl, ctrlEvent, targetFunct ):
        if ctrl is not None:
            self.ctrl           = ctrl
            self.ctrlEvent      = ctrlEvent
            self.targetFunct    = targetFunct

            if self.val is not None:
                self.setVal( self() )

            if self.ctrlEvent is not None:
                self.ctrl.Bind( self.ctrlEvent, self.onPropertyChanged )
        else:
            print("no ctrl given")

    def onPropertyChanged( self, event = None ):
        
        #try:
            #event.Skip()
        #except:
            #print "event cannot be skiped", event

        oldVal = self.val
        self.val = self.getValFromCtrl()

        #print "onPropertyChanged( self, event = None ):", event
        #print self.val
        
        if self.val != oldVal:
            if self.parent is not None:
                self.parent.propertyChanged()

            if self.targetFunct is not None:
                self.targetFunct()
        
    def setVal( self, val ):
        '''
            Convert a value INTO the appropriate ctrl-setting
        '''
#        print "set value for: ", self.name, self.ctrl, val
        if val is None:
            sys.stderr.write( "internal: No value given for setter: " + self.name + self.ctrl )
            return

        if self.ctrl is not None:
            if isinstance( self.ctrl, wx.SpinCtrl ):
                self.ctrl.SetValue( int( val ) )

            elif isinstance( self.ctrl, wx.Choice ):
                self.ctrl.SetSelection( self.ctrl.FindString( val ) )

            elif isinstance( self.ctrl, wx.RadioBox ):
                if self.valType == str or self.valType == str:
                    self.ctrl.SetSelection( self.ctrl.FindString( val ) )
                else:
                    self.ctrl.SetSelection( int( self.val ) )

            elif isinstance( self.ctrl, wx.TextCtrl ):
                if not isinstance( val, str ):
                    self.ctrl.SetValue( str( val ) )
                else:
                    self.ctrl.SetValue( val )

            elif isinstance( self.ctrl, wx.FilePickerCtrl ):
                print(self.ctrl, val)
                self.ctrl.SetPath( val )
            
            elif hasattr( self.ctrl, 'SetValue' ):
                if hasattr( self.ctrl, 'IsChecked' ):
                    self.ctrl.SetValue( bool( val ) )
                elif type( self.ctrl.GetValue() ) == self.valType:
                    self.ctrl.SetValue( val )
                else:
                    self.ctrl.SetValue( str( val ) )

            # in every case get the value back from control
            self.val = self.getValFromCtrl()
        else:
            if self.valType is not None:
                if self.valType == tuple:
                    self.val = [int( x__) for x__ in val[1:len(val)-1 ].split(',')]
                else:
                    self.val = self.valType( val )
            else:
                if val == 'True':
                    self.val = True
                elif val == 'False':
                    self.val = False
                else:
                    self.val = val

    def getValFromCtrl( self ):
        '''
            Convert the value FROM the appropriate ctrl-setting
        '''

        #print "getValFromCtrl", self.ctrl

        if isinstance( self.ctrl, wx.SpinCtrl ):
            return int( self.ctrl.GetValue( ) )

        elif isinstance( self.ctrl, wx.RadioBox ):
            if self.valType == str or self.valType == str:
                return self.ctrl.GetString( self.ctrl.GetSelection( ) )
            else:
                return self.ctrl.GetSelection( )

        elif isinstance( self.ctrl, wx.Choice ):
            return self.ctrl.GetStringSelection( )

        elif isinstance( self.ctrl, wx.TextCtrl ):
            if self.valType == str or self.valType == str or self.valType == None:
                return self.ctrl.GetValue( )
            else:
                if len( self.ctrl.GetValue( ) ) == 0:
                    return 0
                else:
                    return self.valType( self.ctrl.GetValue( ) )

        elif isinstance( self.ctrl, wx.FilePickerCtrl ):
            print(self.ctrl, self.ctrl.GetPath())
            return self.ctrl.GetPath()

        elif hasattr( self.ctrl, 'IsChecked' ):
            return self.ctrl.IsChecked();
        elif hasattr( self.ctrl, 'GetValue' ):
            return self.ctrl.GetValue()

        print("Cannot convert Value from control. Pls fix me", self.ctrl)
    
# END class Property
                            
class ManagedProperties:
    """What is this??"""
    def __init__( self ):
        self.properties = dict()
        self.piCaption  = "None"
        self.debug_     = True
            
    def appendProperty( self, name, default = None, valType = None
                        , ctrl = None, ctrlEvent = None, targetFunct = None ):
        """"""
        prop = None
        
        if name in self.properties:
            prop = self.properties[ name ]
        else:
            prop = Property( name, default, valType )
            self.properties[ name ] = prop
            prop.parent = self

        if ctrl is not None:
            prop.setCtrl( ctrl, ctrlEvent, targetFunct )

            if default is None:
                prop.val = prop.getValFromCtrl()

        return prop

    def writePropertiesToFile( self, fi ):
        fi.write( "[" + self.piCaption + "]\n" )

        for k, v in self.properties.items():
            print("write:", fi, k, v, v())
            fi.write( k + "=" + str( v() ) + "\n" )

    def setProperties( self, props ):
        #print "setProperties:", self.piCaption, props
        for k, v in props.items():
            if k in self.properties:
                print(k, v, self.properties[ k ])
                self.properties[ k ].setVal( v )
            else:
                sys.stderr.write( "unknown property for panel "  + self.piCaption + " : " + k + "\n" )

    def propertyChanged( self ):
        ''' Abstract: will be called if some property value changed. to be overridden by child classes '''
        pass

# END class ManagedProperties
        
class AppResource( ManagedProperties ):
    """What is this??"""
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        ManagedProperties.__init__( self )
        
        self.xrc                 = None
        self.parent              = parent

        self.name_                   = "No Name"
        self.rendererSlot_           = rendererSlot
        self.propertyInspectorSlot_  = propertyInspectorSlot

        self.propertyPanel_     = None
        self.rendererPanel_     = None
        self.renderer_          = None

        self.treeItem           = None
        self.active             = False
        self.parentResource     = None

        self.dependencyList     = []
        self.listenerList       = []

        self.dataSources_       = dict()

        self.subPanels          = []

        self.mainMenus          = dict()
        
        # store all local menus here
        self.menuSections       = []

        self.dependencyChanged_ = True

    def setName( self, name ) : self.name_ = name;
    
    def getName( self ) : return self.name_;
    
    def getResourceTree( self ): return self.parent.resourceTree

    def loadXRC( self, xrcfile, localfile ):
        self.xrc = xrc.EmptyXmlResource()

        if hasattr( sys, "frozen"):
            globPath = os.path.dirname( sys.argv[ 0 ] )
        else:
            globPath = os.path.dirname( localfile )

        mainXRCFile=os.path.join( globPath, xrcfile )

        if not os.path.exists( mainXRCFile ):
            raise IOError('Could not find xrc file "%s"; dying'%mainXRCFile)

        self.xrc.Load( mainXRCFile )

    def destroy( self ):
        self.activate( False )
        self.activateApplication( False )

    def setSource( self, name, function ):
        """Add a functor to source data into the data dictionary."""
        self.dataSources_[ name ] = function

    def getSource( self, name ):
        """Get source data from the functor the data dictionary."""
        if name in self.dataSources_:
            return self.dataSources_[ name ]()
        else:
            raise Exception( 'no data source defined for: ' + str( name ) )

    def getDependency( self, classname ):
        for c in self.dependencyList:
            if isinstance( c, classname ):
                return c

    def addDependency( self, c ):
        if c not in self.dependencyList:
            self.dependencyList.append( c )

        if self not in c.listenerList:
            c.listenerList.append( self )

    def dependencyChanged( self ):
        print("dependencyChanged( self )", self)
        self.dependencyChanged_ = True

    def notifyListenerDependencyChanged_( self, force = False ):
        for l in self.listenerList:
            l.dependencyChanged()
            #if force:
                #l.process( False )
                #return
            #else:
                #if hasattr( l, 'IsShownOnScreen' ):
                    #if l.IsShownOnScreen():
                        #l.process( False )
                        #return

    def checkDependencies( self, recursive = True ):
        print("CheckDependencies:", self)
        print(self.listenerList)
        if recursive:
            for d in self.dependencyList:
                d.checkDependencies( )

        if self.processData():
            self.notifyListenerDependencyChanged_()

        self.dependencyChanged_ = False

    def process( self, recursive = True ):
        print("process:", self)
        self.checkDependencies( recursive )

        # update resources
        if self.parentResource:
            for s in self.parentResource.subPanels:
                if s is not self:
                    print(self, s, s.dependencyChanged_, s.IsShownOnScreen())
                    if s.dependencyChanged_ and s.IsShownOnScreen():
                        s.process()

            if self.parentResource.dependencyChanged_ and not hasattr( self.parentResource, 'IsShownOnScreen' ):
                print("parent process forced:", self.parentResource)
                self.parentResource.process( recursive = False )


    def processData( self ):
        """
        To be overridden by child.

        Return True if something was done, that the listener should
        know, else return False
        """
        return False;

    def getPropertyPanel( self, slot ):
        """What is this?"""
        if self.propertyPanel_ is None and slot is not None:
            if hasattr( self, 'createPropertyPanel' ):
                self.propertyPanel_ = self.createPropertyPanel( slot )
                slot.GetSizer().Add( self.propertyPanel_, 1, wx.EXPAND )
                slot.GetSizer().Layout()

        return self.propertyPanel_

    def getRenderer( self, slot = None):
        """What is this?"""
        
        if not self.renderer_ and slot is not None:
            self.renderer_ = aui.AuiNotebook( slot
                                            , style = aui.AUI_NB_TOP | aui.AUI_NB_TAB_SPLIT | aui.AUI_NB_TAB_MOVE | aui.AUI_NB_SCROLL_BUTTONS )

            if hasattr( self, 'createRendererPanel' ):
                self.rendererPanel_ = self.createRendererPanel( slot )
                self.renderer_.AddPage( self.rendererPanel_, self.rendererPanel_.GetName(), False)

            self.renderer_.Bind( aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.onRendererTabSwitch )
        else:
            print("def getRenderer( self, slot = None):")

        print(self.renderer_, slot)
        return self.renderer_

    def getRendererPanel( self ) : return self.rendererPanel_

    def onRendererTabSwitch( self, event ):
        ''
        '  What happens if a tab from the renderview is switched '
        ''
        newTab = self.renderer_.GetPage( event.GetSelection() )

        if newTab == self.rendererPanel_:
            newTab = self

        #print "onRendererTabSwitch", newTab
        
        # change the selected tree item 
        self.parent.resourceTree.selectItem( newTab )

        event.Skip()
    # def onRendererTabSwitch( ... )
    
    def activatePropertyPanel( self, active ):
        """What is this?"""
        p = self.getPropertyPanel( self.propertyInspectorSlot_ )
        if p is not None:
            if active:
                p.Show()
            else:
                p.Hide()

    def activateToolBar_( self, tb, active, pos = 5):
        """What is this?"""
        #print "activateToolBar( self, active ):", self, active
        if active:
            name = tb.GetName()
            if len( name ) == 0:
                name = "Toolbar"
        #print "add pane", name
            self.parent.auiMgr.AddPane( tb, aui.AuiPaneInfo()
                                            .Name( name ).Caption( name )
                                            .ToolbarPane().Top().Row(0).Position(pos)
                                            .LeftDockable(True).RightDockable(True)
            #self.parent.auiMgr.InsertPane( tb, pos
                                            )
            tb.Show()
        else:
            #print "del pane", tb
            tb.Hide()
            self.parent.auiMgr.DetachPane( tb )

        self.parent.auiMgr.Update()

    def activateApplication( self, active ):
        ''
        ' What is this? '
        ''
        print("activateApplication( self, active ):", self, active)
        
        if hasattr( self, "getRenderer" ):
            r = self.getRenderer( self.rendererSlot_ )
            if active:
                print("activate parent: ", r)
                self.rendererSlot_.GetSizer().Add( r, 1, wx.EXPAND, 0  )
                self.rendererSlot_.GetSizer().Layout()
                r.Show()

                #if hasattr( getApplicationToolBar
            else:
                #print "deactivate parent: ", r
                r.Hide()
                self.rendererSlot_.GetSizer().Detach( r )
                self.rendererSlot_.GetSizer().Layout()
                    #if self.parentResource:
            #self.parentResource.activatePropertyPanel( False )
            #if self.parentResource.active != active:
                #return
        
        if hasattr( self, 'getApplicationToolBar' ):
            self.activateToolBar_( self.getApplicationToolBar( self.parent ), active, pos = 1 )

        for m in list(self.mainMenus.keys()):
            pos = self.parent.GetMenuBar().FindMenu( m.GetTitle() )
            #print "check: ", m.GetTitle(), pos
            if active:
                if pos is wx.NOT_FOUND:
                    #print "activate"
                    self.parent.GetMenuBar().Insert( 2, m, self.mainMenus[ m ] )
                    m.SetTitle( self.mainMenus[ m ].replace( '&', ''))
            else:
                if pos is not wx.NOT_FOUND:
                    #print "deactivate"
                    self.parent.GetMenuBar().Remove( pos )

    def activate( self, active ):
        ''
        ' What is this? '
        ''
        print("activate( self, active ):", self, self.active, active)
        
        if self.active == active:
            return

        self.active = active

        #if activates rise appropriate render tab, if available
        if self.active:
            nb = None

            if isinstance( self.rendererSlot_, aui.AuiNotebook ):
                print("# self is child renderpanel (AppResourceWx*)")
                nb = self.rendererSlot_
            elif isinstance( self.getRenderer( ), aui.AuiNotebook ):
                print("# self is parent (AppResource)")
                nb = self.getRenderer( )

            if nb is not None:
                parentRendererPanel = nb.GetPageIndex( self.getRendererPanel() )
        #print "parentRendererPanel", parentRendererPanel
                if nb.GetSelection() != parentRendererPanel and parentRendererPanel > -1:
                    nb.SetSelection( parentRendererPanel )

        if hasattr( self, 'getToolBar' ):
            self.activateToolBar_( self.getToolBar( self.parent ), active )
            
        self.activatePropertyPanel( active )
        
        # show/hide Local menusection
        for sec in self.menuSections:
            sec.activate( active )
        
        if self.active and self.dependencyChanged_:
            print("process due to activate", self)
            self.process()
            
        self.parent.resourceTree.selectItem( self )
        

        self.parent.auiMgr.Update()

    def createSubPanel( self, classname, name = None ):
        """What is this?"""
        r = self.getRenderer( self.rendererSlot_ )
        
        if isinstance( r , aui.AuiNotebook ):
            panel = classname( self.parent, r, self.propertyInspectorSlot_ );
            if panel.xrc is None:
                panel.xrc = self.xrc

            if name is None:
                name = panel.getName()

            r.AddPage( panel, name, False );

            self.parent.resourceTree.addItem( panel, name, self.treeItem )
            panel.parentResource = self
            self.subPanels.append( panel )
            return panel
        else:
            raise Exception( "createSubPanel only defined for application with notebook renderer" )

    def propertyChanged( self ):
        """What is this."""
        if self.parentResource is not None:
            self.parentResource.propertyChanged()

    def createPropertyInspectorNoteBookPanel( self, parent, panelName, title = None ):
        '''
            Create and return the Property-Notebook and create the first tab from panelName
        '''
        if title is not None:
            self.piCaption = title

        propertyNoteBook = aui.AuiNotebook( parent,
                                            style = aui.AUI_NB_TOP | aui.AUI_NB_TAB_SPLIT |
                                                    aui.AUI_NB_TAB_MOVE | aui.AUI_NB_SCROLL_BUTTONS )

        scrolled = wx.lib.scrolledpanel.ScrolledPanel( propertyNoteBook
                                        , -1, style = wx.TAB_TRAVERSAL, name="PropertySlot" + panelName ) 
        
        #print self.xrc
        #print self.parent.xrc

        tab = None
        if self.xrc is not None:
        #    print "self.xrc "
            tab = self.xrc.LoadPanel( scrolled, panelName )
        else:
         #   print "self.parent.xrc "
            tab = self.parent.xrc.LoadPanel( scrolled, panelName )
        
        sizer = wx.BoxSizer()
        sizer.Add( tab, 1 )
        sizer.Layout()
        
        scrolled.SetSizer( sizer )
        scrolled.SetAutoLayout(1)
        scrolled.SetupScrolling()

        scrolled.Layout()
        
        if scrolled is not None:
            propertyNoteBook.AddPage( scrolled, self.piCaption, True )
            
        return propertyNoteBook
        
    def setStatusMessage( self, msg ):
        """Helper function for manageing status messages."""
        self.parent.statusBar.setStatusMessage( msg )
    #def setStatusMessage( ... )

#    
#START Menu related stuff    
#
    def findMainMenu( self, name ):
        """
        Find the main menu entry with a given name.

        If none exist create one.
        """
        self.mainMenus          = dict()
        pos = self.parent.GetMenuBar().FindMenu( name )
        if pos > -1:
            return self.parent.GetMenuBar().GetMenu( pos )
        else:
            return self.createMainMenu( name )
    # def findMainMenu( ... ):
    
    def createMainMenu( self, name, pos = 3 ):
        """Create new main menu entry."""
        menu = wx.Menu( )
        self.parent.GetMenuBar().Insert( pos, menu, name )
        menu.SetTitle( name.replace( '&', '' ) )
        self.mainMenus[ menu ] = name
        
        return menu
        
    def createLocalMenuSection( self, menu ):
        sec = MenuSection( menu )
        self.menuSections.append( sec )
        return sec

    def createMenuItem( self, menu, name = "", help = "auto", function = None, bitmap = None ):
        """
        What is this?

        If help is set to 'auto', the docstring of function is used. please ensure this docstring is a one-liner
        """
        if help == 'auto':
            help = function.__doc__
        item = wx.MenuItem( menu, wx.NewId(), name, help )
        
        if bitmap is not None:
            item.SetBitmap( bitmap )
            
        menu.AppendItem( item )
        
        self.parent.Bind( wx.EVT_MENU, function, id = item.GetId() )
        return item
        
    def getExportFileNameWithDialog( self, defaultFile, wildcard, message = "Choose file"
                                   , defaultDir = os.getcwd()
                                   , style=None):
        if not style:
            style = wx.SAVE|wx.OVERWRITE_PROMPT|wx.CHANGE_DIR
            
        dlg = wx.FileDialog( self.parent
                            , message = message
                            , defaultDir = defaultDir
                            , defaultFile = defaultFile
                            , wildcard = wildcard
                            , style = style )
                            
        #dlg.SetFilterIndex(filter_index)
        if dlg.ShowModal() == wx.ID_OK:
            paths       = dlg.GetPaths()
            dirname     = dlg.GetDirectory()
            filenames    = dlg.GetFilenames()
            
            print("getExportFileNameWithDialog", paths, dirname, filenames)
            
            if style & wx.FD_MULTIPLE:
                return [unicodeToAscii( f_ ) for f_ in paths]
            else:
                return unicodeToAscii( paths[ 0 ] )

        return None

#    
#END Menu related stuff    
#

# END class AppResourceWx

class MenuSection:
    """This is a helper class to manage local menu groups, that can be
    shown/hidden by calling activate."""
    def __init__( self, menu ):
        self.menu = menu
        self.items = []
        self.names = []
        self.targetFunctions = []
        
    def addItem( self, name = "", help = "auto", function = None, bitmap = None ):
        """"""
        if help == 'auto':
            help = function.__doc__
            
        item = wx.MenuItem( self.menu, wx.NewId(), name, help )
        self.items.append( item )
        self.names.append( name )
        self.targetFunctions.append( function )
        
        if bitmap is not None:
            item.SetBitmap( bitmap )

        self.menu.Bind( wx.EVT_MENU, function, id = item.GetId() )
    #def addItem( ... )
            
    def isItemInMenu( self, item, menu ):
        #print "############################################################"
        #print item
        #print item.GetLabel()
        #print "############################################################"
        for i in range( menu.GetMenuItemCount() ):
            #testItem = menu.FindItemByPosition( i )
            #print "test:", item.GetLabel(), " vs. ", testItem.GetLabel()
            
        
            #print "item:", item, id( item ), item.Id
            #print "find:", menu.FindItemByPosition( i ), id(menu.FindItemByPosition( i )), menu.FindItemByPosition( i ).Id
            #print item is menu.FindItemByPosition( i )
            #print item == menu.FindItemByPosition( i )
            #print item.Id == menu.FindItemByPosition( i ).Id
            if item.Id == menu.FindItemByPosition( i ).Id: return True
                
        return False

    def activate( self, activate ):
        if activate:
            # add two seps (leeding and trailing); remove the trailing later, while .AppendItem segfaults (see below)
            self.menu.AppendSeparator( )
            self.menu.AppendSeparator( )
            
            # dämlicher Check weil wx keine gescheite funktion anbietet
            for item in self.items:
                if not self.isItemInMenu( item, self.menu ):
                    #print "Activate: 1 ", item, self.menu.GetMenuItemCount()
                    self.menu.InsertItem( self.menu.GetMenuItemCount()-1, item )
                    
                    # append function segfaults here within this testing context, so we need a evil hack
                    #print "Activate: 1 ", item, self.menu.GetMenuItemCount()
                    #self.menu.AppendItem( item ) # fails
                    #print "Activate: 2 ", item, self.menu.GetMenuItemCount()
                    #self.menu.RemoveItem( item )
                    #print "Activate: 3 ", item, self.menu.GetMenuItemCount()
                    #self.menu.AppendItem( item ) # fails
                    #print "Activate: 4 ", item, self.menu.GetMenuItemCount()

            #remove trailing here
            self.menu.RemoveItem( self.menu.FindItemByPosition( self.menu.GetMenuItemCount() -1 ) )
        else:
            for item in self.items:
                if self.isItemInMenu( item, self.menu ):
                    #print "Activate: remove", item
                    self.menu.RemoveItem( item )
                    #print "Activate: removed", item
                    
            #remove leeding sep here
            self.menu.RemoveItem( self.menu.FindItemByPosition( self.menu.GetMenuItemCount() -1 ) )

#end class MenuSection

