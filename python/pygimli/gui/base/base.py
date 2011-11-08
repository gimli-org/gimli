# -*- coding: utf-8 -*-

import sys, os
import wx
import wx.xrc as xrc

try:
    from agw import aui
    from agw.aui import aui_switcherdialog as ASD
except ImportError: # if it's not there locally, try the wxPython lib.
    import wx.lib.agw.aui as aui
    from wx.lib.agw.aui import aui_switcherdialog as ASD

from pygimli.utils import unicodeToAscii

class Property:
    '''
        What is this??
    '''
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
            print "no ctrl given"

    def onPropertyChanged( self, event = None ):
        try:
            event.Skip()
        except:
            pass
        
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
                if self.valType == str or self.valType == unicode:
                    self.ctrl.SetSelection( self.ctrl.FindString( val ) )
                else:
                    self.ctrl.SetSelection( int( self.val ) )

            elif isinstance( self.ctrl, wx.TextCtrl ):
                if not isinstance( val, basestring ):
                    self.ctrl.SetValue( str( val ) )
                else:
                    self.ctrl.SetValue( val )

            elif isinstance( self.ctrl, wx.FilePickerCtrl ):
                print self.ctrl, val
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
                    self.val = map( lambda x__: int( x__), val[1:len(val)-1 ].split(',') )
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
            if self.valType == str or self.valType == unicode:
                return self.ctrl.GetString( self.ctrl.GetSelection( ) )
            else:
                return self.ctrl.GetSelection( )

        elif isinstance( self.ctrl, wx.Choice ):
            return self.ctrl.GetStringSelection( )

        elif isinstance( self.ctrl, wx.TextCtrl ):
            if self.valType == str or self.valType == unicode or self.valType == None:
                return self.ctrl.GetValue( )
            else:
                if len( self.ctrl.GetValue( ) ) == 0:
                    return 0
                else:
                    return self.valType( self.ctrl.GetValue( ) )

        elif isinstance( self.ctrl, wx.FilePickerCtrl ):
            print self.ctrl, self.ctrl.GetPath()
            return self.ctrl.GetPath()

        elif hasattr( self.ctrl, 'IsChecked' ):
            return self.ctrl.IsChecked();
        elif hasattr( self.ctrl, 'GetValue' ):
            return self.ctrl.GetValue()

        print "Cannot convert Value from control. Pls fix me", self.ctrl
    
# END class Property
                            
class ManagedProperties:
    '''
        What is this??
    '''
    def __init__( self ):
        self.properties = dict()
        self.piCaption  = "None"
        self.debug_     = True
            
    def appendProperty( self, name, default = None, valType = None
                        , ctrl = None, ctrlEvent = None, targetFunct = None ):
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

        for k, v in self.properties.iteritems():
            print "write:", fi, k, v, v()
            fi.write( k + "=" + unicode( v() ) + "\n" )

    def setProperties( self, props ):
        #print "setProperties:", self.piCaption, props
        for k, v in props.iteritems():
            if k in self.properties:
                print k, v, self.properties[ k ]
                self.properties[ k ].setVal( v )
            else:
                sys.stderr.write( "unknown property for panel "  + self.piCaption + " : " + k + "\n" )

    def propertyChanged( self ):
        ''' Abstract: will be called if some property value changed. to be overridden by child classes '''
        pass

# END class ManagedProperties
        
class AppResource( ManagedProperties ):
    '''
        What is this??
    '''
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
        self.dataSources_[ name ] = function

    def getSource( self, name ):
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
        print "dependencyChanged( self )", self
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
        print "CheckDependencies:", self
        print self.listenerList
        if recursive:
            for d in self.dependencyList:
                d.checkDependencies( )

        if self.processData():
            self.notifyListenerDependencyChanged_()

        self.dependencyChanged_ = False

    def process( self, recursive = True ):
        print "process:", self
        self.checkDependencies( recursive )

        # update resources
        if self.parentResource:
            for s in self.parentResource.subPanels:
                if s is not self:
                    print self, s, s.dependencyChanged_, s.IsShownOnScreen()
                    if s.dependencyChanged_ and s.IsShownOnScreen():
                        s.process()

            if self.parentResource.dependencyChanged_ and not hasattr( self.parentResource, 'IsShownOnScreen' ):
                print "parent process forced:", self.parentResource
                self.parentResource.process( recursive = False )


    def processData( self ):
        '''
            To be overridden by child. Return True if something was done, that the listener should know, else return False
        '''
        return False;

    def createMainMenu( self, name, pos = 3 ):
        '''
            What is this?
        '''
        menu = wx.Menu( )
        print "create createMainMenu", name 
        #inserted by activateApplication
        self.parent.GetMenuBar().Insert( pos, menu, name )
        menu.SetTitle( name.replace( '&', '' ) )
        self.mainMenus[ menu ] = name
        return menu

    def createMenuItem( self, menu, name = "", help = "", function = None, bitmap = None ):
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

    def getPropertyPanel( self, slot ):
        '''
            What is this?
        '''
        if self.propertyPanel_ is None and slot is not None:
            if hasattr( self, 'createPropertyPanel' ):
                self.propertyPanel_ = self.createPropertyPanel( slot )
                slot.GetSizer().Add( self.propertyPanel_, 1, wx.EXPAND )
                slot.GetSizer().Layout()

        return self.propertyPanel_

    def getRenderer( self, slot = None):
        '''
            What is this?
        '''
        
        if not self.renderer_ and slot is not None:
            self.renderer_ = aui.AuiNotebook( slot
                                            , style = aui.AUI_NB_TOP | aui.AUI_NB_TAB_SPLIT | aui.AUI_NB_TAB_MOVE | aui.AUI_NB_SCROLL_BUTTONS )

            if hasattr( self, 'createRendererPanel' ):
                self.rendererPanel_ = self.createRendererPanel( slot )
                self.renderer_.AddPage( self.rendererPanel_, self.rendererPanel_.GetName(), False)

            self.renderer_.Bind( aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.onRendererTabSwitch )
        else:
            print "def getRenderer( self, slot = None):"

        print self.renderer_, slot
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
        '''
            What is this?
        '''
        p = self.getPropertyPanel( self.propertyInspectorSlot_ )
        if p is not None:
            if active:
                p.Show()
            else:
                p.Hide()

    def activateToolBar_( self, tb, active, pos = 5):
        '''
            What is this?
        '''
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
        print "activateApplication( self, active ):", self, active
        
        if hasattr( self, "getRenderer" ):
            r = self.getRenderer( self.rendererSlot_ )
            if active:
                print "activate parent: ", r
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

        for m in self.mainMenus.keys():
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
        print "activate( self, active ):", self, self.active, active
        
        if self.active == active:
            return

        self.active = active

        #if activates rise appropriate render tab, if available
        if self.active:
            nb = None

            if isinstance( self.rendererSlot_, aui.AuiNotebook ):
                print "# self is child renderpanel (AppResourceWx*)"
                nb = self.rendererSlot_
            elif isinstance( self.getRenderer( ), aui.AuiNotebook ):
                print "# self is parent (AppResource)"
                nb = self.getRenderer( )

            if nb is not None:
                parentRendererPanel = nb.GetPageIndex( self.getRendererPanel() )
        #print "parentRendererPanel", parentRendererPanel
                if nb.GetSelection() != parentRendererPanel and parentRendererPanel > -1:
                    nb.SetSelection( parentRendererPanel )

        if hasattr( self, 'getToolBar' ):
            self.activateToolBar_( self.getToolBar( self.parent ), active )
            
        self.activatePropertyPanel( active )
        
        if self.active and self.dependencyChanged_:
            print "process due to activate", self
            self.process()
            
        self.parent.resourceTree.selectItem( self )
        

        self.parent.auiMgr.Update()

    def createSubPanel( self, classname, name = None ):
        '''
            What is this?
        '''
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
        '''
            What is this
        '''
        if self.parentResource is not None:
            self.parentResource.propertyChanged()

    def createPropertyInspectorNoteBookPanel( self, parent, panelName, title = None ):
        '''
            What is this
        '''
        if title is not None:
            self.piCaption = title

        panel = aui.AuiNotebook( parent
                                            , style = aui.AUI_NB_TOP | aui.AUI_NB_TAB_SPLIT |
                                            aui.AUI_NB_TAB_MOVE | aui.AUI_NB_SCROLL_BUTTONS )

        #print self.xrc
        #print self.parent.xrc
        tab = None
        if self.xrc is not None:
        #    print "self.xrc "
            tab = self.xrc.LoadPanel( parent, panelName )
        else:
         #   print "self.parent.xrc "
            tab = self.parent.xrc.LoadPanel( parent, panelName )

        if tab is not None:
            panel.AddPage( tab, self.piCaption, True )
        return panel
# END class AppResource

