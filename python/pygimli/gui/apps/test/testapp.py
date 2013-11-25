# -*- coding: iso-8859-1 -*-

import wx

from pygimli.gui.base import AppResource

try:
    from pygimli.gui.wxmpl import AppResourceWxMPL
except Exception as e:
	print 'Import of wxmpl failed: ', e

try:
    from pygimli.gui.wxvtk import vtk, AppResourceWxVTK
except:
    pass

from pygimli.gui.resources import loadIcon, MakeDisabledBitmap

HAVE_PROPGRID = False

class TestVTKApp( AppResourceWxVTK ):
    """this is a vtk child renderer for the test application see below."""
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        AppResourceWxVTK.__init__( self, parent, rendererSlot, propertyInspectorSlot )
        
        cone = vtk.vtkConeSource()
        cone.SetResolution(8)
    
        coneMapper = vtk.vtkPolyDataMapper()
        coneMapper.SetInput(cone.GetOutput())
    
        coneActor = vtk.vtkActor()
        coneActor.SetMapper(coneMapper)

        textActor = vtk.vtkTextActor()
        textActor.ScaledTextOn()
        textActor.SetDisplayPosition(90, 50)
        textActor.SetInput("This is a text example")
        
        # Set coordinates to match the old vtkScaledTextActor default value
        textActor.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
        textActor.GetPosition2Coordinate().SetValue(0.6, 0.1)

        tprop = textActor.GetTextProperty()
        tprop.SetFontSize(18)
        tprop.SetFontFamilyToArial()
        tprop.SetJustificationToCentered()
        tprop.BoldOff()
        tprop.ItalicOff()
        tprop.ShadowOff()
        tprop.SetColor(0, 0, 0)

        self.renderer.AddActor2D(textActor)
        self.renderer.AddActor(coneActor)

    def createPropertyPanel( self, parent ):
        if not HAVE_PROPGRID:
            return wx.Panel( parent )
        panel = wxPropGridWrapper( parent )
        panel.Append( wxpg.PropertyCategory("viewport") )
        return panel

class TestMPLApp( AppResourceWxMPL ):
    """this is a mpl child renderer for the test application see below."""
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        AppResourceWxMPL.__init__( self, parent, rendererSlot, propertyInspectorSlot )
            
    def createPropertyPanel( self, parent ):
        if not HAVE_PROPGRID:
            return wx.Panel( parent )
    
        panel = wxPropGridWrapper( parent )
        
        panel.Append( wxpg.PropertyCategory("xaxis") )
        self.fillPGWithMPLObject( panel, self.axes.xaxis, self.onXAxisChanged )
        panel.Append( wxpg.PropertyCategory("yaxis") )
        self.fillPGWithMPLObject( panel, self.axes.yaxis, self.onYAxisChanged )
        panel.Append( wxpg.PropertyCategory("axes") )
        self.fillPGWithMPLObject( panel, self.axes.axes, self.onAxesChanged )
        
        return panel
        
    def fillPGWithMPLObject( self, panel, obj, target ):
        for att in dir( obj ):
            try:
                if hasattr( getattr( obj, att), '__call__' ) and  \
                    att.find('set_') == 0 :
                    
                    getter = att.replace( 'set_', 'get_' )
                    if hasattr( getattr( obj, getter ), '__call__' ):
                        retVal = getattr( obj, getter )()
                        
                        #print retVal, type( retVal )
                        if type( retVal ) is type( True ):
                            prop = wxpg.BoolProperty( att, value = retVal )
                        elif type( retVal ) is type( 1.0 ):
                            prop = wxpg.FloatProperty( att, value = retVal )
                        else:
                            prop = wxpg.StringProperty( att, value = str( retVal ) )
                            prop.disable()
                        panel.append( prop, target, toolTip = obj.__getattribute__( att ).__doc__ )
            except:
                pass
        #panel.append( wxpg.StringProperty( "x-label", value = "name" ), self.onXLableChanged, toolTip = "name" )
        
    def onXAxisChanged( self, name, value ):
        getattr( self.axes.xaxis, name )( value )
        self.updateDrawOnIdle()
        
    def onYAxisChanged( self, name, value ):
        getattr( self.axes.xaxis, name )( value )
        self.updateDrawOnIdle()
        
    def onAxesChanged( self, name, value ):
        getattr( self.axes, name )( value )
        self.updateDrawOnIdle()

    def onMousePress( self, event ):
        print "mousePress"
        axes = event.inaxes
        if axes is not None:
            #print axes
            print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)
            p = axes.pick(event)
            print "inaxes pick", p
            
        else:
            print "Pos:", event.canvas._lastx, event.canvas._lasty
            xdata, ydata = self.axes.transData.inverted().transform_point((event.canvas._lastx, event.canvas._lasty))
            print "DPos:", xdata, ydata

    def onPick( self, event ):
        print "onPick", event
        
class TestApp( AppResource ):
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        AppResource.__init__( self, parent, rendererSlot, propertyInspectorSlot )
        
        # set the name of this application that appears to the resource tree
        self.setName("Tester")
        
        ## optionaly load content informations for this application, generated by an gui-builder that wrote xrc-files
        ## eg. wxFormBuilder, __file__ is needed to find the xrc-file 
        self.loadXRC( 'testapp.xrc', __file__ )
        
        ## register a menu to the main-MenuBar, using this function ensures that this menu will be shown/hiden respective to the activate status of the application
        self.mbTestMenu = self.createMainMenu( '&Test' )
        
        ## create an item to to the test menu
        self.createMenuItem( self.mbTestMenu
                                    , name = "Test Menu item"
                                    , help = "Test Menu item"
                                    , function = self.menuFunction )
                                    
        ## the application has some properties that can be altered by the property inspector (PI), loaded and saved
        self.titleTextProp = self.appendProperty( "Title", default = 'unknown', valType = str )
                
        # do some stuff after this application has created and registered to the resourcetree
        self.parent.addCommandToOnIdleQueue( self.postCreate )
    
    def postCreate( self ):
        self.createSubPanel( TestMPLApp, "wxMPL" )
        self.createSubPanel( TestVTKApp, "wxVTK" )
        
    def createPropertyPanel( self, parent ):
        """
            Define and return panel that is shown in the property-inspector (PI)
            The panel will post created at the first call
        """

        # create a Notebook for the PI and add the content for the panel with the name piTestApplication defined in testapp.xrc
        panel = self.createPropertyInspectorNoteBookPanel( parent, 'piTestApplication', title = 'TestApplication' )
            
        # define property behaviour
        self.titleTextProp.setCtrl( ctrl = wx.xrc.XRCCTRL( panel,  'TitleTextCtrl' ) # name of the control in xrc
                                        , ctrlEvent = wx.EVT_TEXT                        # the event that should observed
                                        , targetFunct = self.setTitle )                  # the callback when the event is called
        
        return panel
        
    def createRendererPanel( self, parent ):
        """Optionally create a renderer panel for this application beside the
        subpanels The panel will post created at the first call."""
        panel = wx.TextCtrl( parent, -1, "", 
                            wx.DefaultPosition, wx.DefaultSize, wx.TE_MULTILINE|wx.NO_BORDER )
        panel.SetName( "Log or Whatever" )
            #S = scrolled.ScrolledPanel( self.rendererPanel_ )
            #S.set_policy(gtk.POLICY_AUTOMATIC,gtk.POLICY_AUTOMATIC)
            #self.ipython = IPythonView()
            ##V.modify_font(pango.FontDescription(FONT))
            #self.ipyton.set_wrap_mode(gtk.WRAP_CHAR)
            #self.ipyton.show()
            #S.add( self.ipython )

        return panel
    
    def menuFunction( self, event = None ):
        """do nothing here, just to show the technique."""
        pass
    
    def openFile( self, files = None ):
        """do nothing here, just to show the technique for MainOpenFileSuffix
        register."""
        print files
    
    def setTitle( self, title = None ):
        """Set the title, either called manually or by callback."""
        if title is None:
            title = self.titleTextProp()
        else:
            self.titleTextProp.setVal( title )

        self.getResourceTree().SetItemText( self.treeItem, self.getName() + ": " + title )

