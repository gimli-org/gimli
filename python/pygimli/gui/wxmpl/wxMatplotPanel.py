# -*- coding: utf-8 -*-
# taken and modified from http://www.scipy.org/Matplotlib_figure_in_a_wx_panel

import sys
import os

try:
    import wx
    import wx.lib.scrolledpanel as scrolled
except ImportError as e:
    import traceback
    #traceback.print_exc(file=sys.stdout)
    sys.stderr.write("No proper wx installed'.\n")

import matplotlib

#matplotlib.interactive( False )
matplotlib.use( 'WXAgg' )

import matplotlib.cbook as cbook

from pygimli.gui.resources import loadIcon
from pygimli.gui.resources import MakeDisabledBitmap
from pygimli.gui.patches import PatchedAuiToolBar
from pygimli.gui.controls import ListCtrlComboPopupControl
from pygimli.gui.base import AppResource
from pygimli.gui.base import AppResource

from pygimli import Stopwatch

try:
    from agw import aui
    from agw.aui import aui_switcherdialog as ASD
except ImportError: # if it's not there locally, try the wxPython lib.
    import wx.lib.agw.aui as aui
    from wx.lib.agw.aui import aui_switcherdialog as ASD

class MPLRubberBander:
    """
    A class to manage mouse events / rubberbanding of a mpl canvas object.

    Instance must be started and stoped
    """
    def __init__(self, canvas ):

        # canvas object
        self.canvas_ = canvas
        # mouse selection start point
        self.startPoint_ = wx.Point(0,0)
        # mouse selection end point
        self.endPoint_   = wx.Point(0,0)
        # mouse selection cache point
        self.savePoint_  = wx.Point(0,0)

        # flags for left click/ selection
        self.button_ = 0
        self.rectShown_ = False

        self.lastAxes_ = None

        self.mouseMoveID_ = None
        self.mousePressID_ = None
        self.mouseReleaseID_ = None
        # Register event handlers for mouse

        self.zoomCallback_ = None
        self.panCallback_ = None

        self.oldCursor_ = None

    def registerEventHandlers(self):
        """Register event handlers for this object."""
        if self.mouseMoveID_ is None:
            self.mouseMoveID_ = self.canvas_.mpl_connect( 'motion_notify_event', self.onMouseMove )
        if self.mousePressID_ is None:
            self.mousePressID_ = self.canvas_.mpl_connect( 'button_press_event', self.onMousePress )
        if self.mouseReleaseID_ is None:
            self.mouseReleaseID_ = self.canvas_.mpl_connect( 'button_release_event', self.onMouseRelease )

    def releaseEventHandlers(self):
        """Release events."""
        if self.mouseMoveID_ is not None:
            self.mouseMoveID_ = self.canvas_.mpl_disconnect( self.mouseMoveID_ )
        if self.mousePressID_ is not None:
            self.mousePressID_ = self.canvas_.mpl_disconnect( self.mousePressID_ )
        if self.mouseReleaseID_ is not None:
            self.mouseReleaseID_ = self.canvas_.mpl_disconnect( self.mouseReleaseID_ )

    def start( self, zoomCallback, panCallback ):
        self.zoomCallback_ = zoomCallback
        self.panCallback_ = panCallback

        self.oldCursor_ = self.canvas_.GetCursor()
        self.canvas_.SetCursor( wx.StockCursor(wx.CURSOR_CROSS))
        self.registerEventHandlers()

    def stop( self ):
        self.zoomCallback_ = None
        self.panCallback_ = None
        self.releaseEventHandlers()
        self.canvas_.SetCursor( self.oldCursor_ )

    def convertCoords( self, x, y ):
        """convert event.coordinates into global widget coordinates."""
        #print x, y, self.canvas_._lastx, self.canvas_._lasty, self.canvas_.GetSize()

        pos = ( min( max( 0, x ), self.canvas_.GetSize()[0] )
              , min( max( 0, self.canvas_.GetSize()[1] - y ), self.canvas_.GetSize()[1] ) )
        #print pos
        return pos

    def onMousePress( self, event ):
        # determine if mouse start within an axes
        if event.inaxes is not None:
            self.lastAxes_ = event.inaxes

        self.button_ = event.button

        if self.button_ == 2:
            self.canvas_.SetCursor( wx.StockCursor(wx.CURSOR_HAND) )
            self.canvas_.CaptureMouse()
            self.canvas_.ReleaseMouse()

        # Left mouse button down, change cursor to
        # something else to denote event capture
        pos = self.convertCoords( event.x, event.y )
        self.startPoint_ = wx.Point( pos[ 0 ], pos[ 1 ] )

        # invalidate current canvas
        self.canvas_.Refresh()

        # cache current position
        self.savePoint_ = self.startPoint_
        self.lastAxes_ = None

    def onMouseRelease( self, event ):

        if event.inaxes is not None:
            self.lastAxes_ = event.inaxes

        selection = self.getCurrentSelection()
        self.clearCurrentSelection()

        if self.button_ == 2:
            self.canvas_.SetCursor( wx.StockCursor(wx.CURSOR_CROSS))
            self.canvas_.CaptureMouse()
            self.canvas_.ReleaseMouse()

        if self.button_ == 3:
            # set selection to whole window for zoom out
            selection = self.getCurrentSelection()

        #selected area to small
        if abs( selection[0]-selection[2] ) < 5 or abs( selection[1]-selection[3] ) < 5:
            self.button_ = 0 # end of clicking
            return

        if self.lastAxes_ is not None:
            leX, leY = self.lastAxes_.transData.inverted().transform_point( self.convertCoords( selection[0], selection[1] ) )
            riX, riY = self.lastAxes_.transData.inverted().transform_point( self.convertCoords( selection[2], selection[3] ) )

            if self.button_ == 1 or self.button_ == 3:

                if self.zoomCallback_ is not None:
                    self.zoomCallback_( ( selection, ( leX,leY,riX,riY), self.lastAxes_ ) )
        else:
            print("MPLRubberBander::onMouseRelease no axes found")

        self.button_ = 0 # end of clicking

    def onMouseMove( self, event ):
        if event is not None:

            # determine if mouse moves through an axes
            if event.inaxes is not None:
                self.lastAxes_ = event.inaxes

            pos = self.convertCoords( event.x, event.y )
            self.endPoint_ = wx.Point( pos[ 0 ], pos[ 1 ] )

            if self.button_ == 1:
                dc = self.getDC()
                # reset dc bounding box
                dc.ResetBoundingBox()
                dc.BeginDrawing()
                w = (self.savePoint_.x - self.startPoint_.x)
                h = (self.savePoint_.y - self.startPoint_.y)

                # To erase previous rectangle
                dc.DrawRectangle(self.startPoint_.x, self.startPoint_.y, w, h)
                # store last hitted axes
                self.savePoint_ = self.endPoint_ # cache current endpoint

                w = (self.endPoint_.x - self.startPoint_.x)
                h = (self.endPoint_.y - self.startPoint_.y)

                # Set clipping region to rectangle corners
                dc.SetClippingRegion(self.startPoint_.x, self.startPoint_.y, w,h)
                dc.DrawRectangle(self.startPoint_.x, self.startPoint_.y, w, h)
                self.rectShown_ = True;
                dc.EndDrawing()

            elif self.button_ == 2:
                x0, y0 = self.lastAxes_.transData.inverted().transform_point( ( 0.0, 0.0 ) )
                x1, y1 = self.lastAxes_.transData.inverted().transform_point( self.endPoint_ - self.startPoint_ )

                if self.panCallback_ is not None:
                    self.panCallback_( ( ( ( x1 - x0 ), ( y0 - y1 ) ), self.lastAxes_ ) )

                self.startPoint_ = self.endPoint_

    def getDC(self):
        """define rectangle drawing properties."""
        # get device context of canvas
        dc= wx.ClientDC(self.canvas_)

        # Set logical function to XOR for rubberbanding
        dc.SetLogicalFunction(wx.XOR)

        # Set dc brush and pen
        # Here I set brush and pen to white and grey respectively
        # You can set it to your own choices

        # The brush setting is not really needed since we
        # dont do any filling of the dc. It is set just for
        # the sake of completion.

        #wbrush = wx.Brush( wx.Colour(255,255,255), wx.TRANSPARENT)
        wbrush = wx.Brush( wx.Colour(30,30,30) )
        wpen = wx.Pen( wx.Colour(200, 200, 200), 1, wx.SOLID)
        dc.SetBrush(wbrush)
        dc.SetPen(wpen)
        return dc

    def getCurrentSelection(self):
        """Return the current selected rectangle."""
        # if there is no selection, selection defaults to
        # current viewport
        left = wx.Point(0,0)
        right = wx.Point(0,0)

        # user dragged mouse to right
        if self.endPoint_.y > self.startPoint_.y:
            right = self.endPoint_
            left = self.startPoint_
        # user dragged mouse to left
        elif self.endPoint_.y < self.startPoint_.y:
            right = self.startPoint_
            left = self.endPoint_

        return (left.x, left.y, right.x, right.y)

    def clearCurrentSelection(self):
        """Clear the current selected rectangle."""
        dc = self.getDC()
        box = self.getCurrentSelection()
        w = box[2] - box[0]
        h = box[3] - box[1]
        dc.SetClippingRegion( box[0], box[1], w, h )

        if self.rectShown_:
            dc.DrawRectangle( box[0], box[1], w, h )
            self.rectShown_ = False;

        # reset selection to canvas size
        self.resetSelection()

    def resetSelection(self):
        """
            Resets the mouse selection to entire canvas.
        """
        self.startPoint_ = wx.Point(0,0)
        sz = self.canvas_.GetSize()
        w,h = sz.GetWidth(), sz.GetHeight()
        self.endPoint_ = wx.Point(w,h)

#END class MPLRubberBander:

#NavigationToolbar2,
#class wxAUIMatplotPanelToolbar( PatchedAuiToolBar ):
class wxAUIMatplotPanelToolbar( aui.AuiToolBar ):
    CUSTOM_PAN_LEFT_ID          = wx.NewId()
    CUSTOM_PAN_RIGHT_ID         = wx.NewId()
    CUSTOM_HOME_ID              = wx.NewId()
    CUSTOM_SAVE_IMAGE_ID        = wx.NewId()
    CUSTOM_INC_CANVAS_WIDTH_ID  = wx.NewId()
    CUSTOM_DEC_CANVAS_WIDTH_ID  = wx.NewId()
    CUSTOM_FIT_CANVAS_ID        = wx.NewId()
    CUSTOM_ASPECT_AUTO_ID       = wx.NewId()
    CUSTOM_CANVAS_WIDTH_COMBO_ID = wx.NewId()
    CUSTOM_PAGE_ZOOM_ID         = wx.NewId()

    canvasWidthScales =      [  1,    2,    4,    8,    16,    1]
    canvasWidthScalesNames = [ '1x         ',
                               '2x         ',
                               '4x         ',
                               '8x         ',
                               '16x       ',
                               'Fit width' ]
    #_absolute_min_size = wx.Size( 100, 22 )
    #_style = AUI_TB_VERTICAL

    def __init__(self, panel, canvas, auiTarget = None):
        self.panel = panel
        
        if not auiTarget:
            auiTarget = panel

#        PatchedAuiToolBar.__init__( self, auiTarget, -1, wx.DefaultPosition, wx.DefaultSize,
#                                    aui.AUI_TB_DEFAULT_STYLE | aui.AUI_TB_PLAIN_BACKGROUND
#                                        #| aui.AUI_TB_TEXT
#                                    )
        aui.AuiToolBar.__init__( self, auiTarget, -1, wx.DefaultPosition, wx.DefaultSize,
                                    aui.AUI_TB_DEFAULT_STYLE | aui.AUI_TB_PLAIN_BACKGROUND
                                        #| aui.AUI_TB_TEXT
                                    )
        self.SetName( "Renderer toolbar" )

        #NavigationToolbar2.__init__(self, canvas)
        self.canvas     = canvas
        self._idle      = True
        self.statbar    = None

        self.SetToolBitmapSize(wx.Size(22,22))

        #self.AddSimpleTool( 'Reset view', self.CUSTOM_HOME_ID, loadIcon( "go-home.png" )
                            #, 'Reset original view' )

        #self.AddCheckTool( self._NTB2_PAN, loadIcon( "my-trans.png" )
                            #, shortHelp='Pan'
                            #, longHelp='Pan with left, zoom with right mouse button')

        #self.AddSimpleTool( self.CUSTOM_PAN_LEFT_ID, loadIcon( "go-previous.png" )
                           #, 'Pan to the left'
                           #, 'Pan view within axes to the left')

        #self.AddSimpleTool( self.CUSTOM_PAN_RIGHT_ID, loadIcon( "go-next.png" )
                           #, 'Pan to the right'
                           #, 'Pan view within axes to the right')

        #self.AddSeparator()

        dec = self.AddSimpleTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, 'Decrease canvas'
                                , loadIcon( "zoom-dec-width.png" )
                                , 'Decrease canvas width' )
        dec.SetDisabledBitmap( MakeDisabledBitmap( loadIcon( "zoom-dec-width.png" ) ) )


        self.canvasComboBox = ListCtrlComboPopupControl( self
                                                        , style=wx.CB_READONLY |wx.CB_DROPDOWN#|wx.NO_BORDER
                                                        , size=(87,-1)
                                                        )
        self.canvasComboBox.popup.dismissCallback = self.onCanvasWidthCombo
        self.canvasComboBox.addItems( self.canvasWidthScalesNames )
        self.canvasComboBox.select( 0 )
        self.canvasComboBox.SetInitialSize( ( 87, 22 ) )
        self.canvasComboBox.SetToolTip( wx.ToolTip( 'Select canvas width zoom' ) )
        self.AddControl( self.canvasComboBox, 'canvas width' )


        inc = self.AddSimpleTool( self.CUSTOM_INC_CANVAS_WIDTH_ID, 'Increase canvas'
                                ,  loadIcon( "zoom-inc-width.png" )
                                , 'Increase canvas width' )
        inc.SetDisabledBitmap( MakeDisabledBitmap(  loadIcon( "zoom-inc-width.png" ) ) )
        self.EnableTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, False )

        self.AddSimpleTool( self.CUSTOM_FIT_CANVAS_ID, 'Fit canvas', loadIcon( "zoom-fit-best.png" )
                           , 'Fit canvas to image size' )

        self.aspectTool = self.AddSimpleTool( self.CUSTOM_ASPECT_AUTO_ID, 'Auto aspect', loadIcon( "transform-scale.png" )
                           , 'Auto aspect', wx.ITEM_CHECK)

        #if ( self.panel.figure.get_axes()[ 0 ].get_aspect() == 'auto' ):
            #self.ToggleTool( self.CUSTOM_ASPECT_AUTO_ID, True )


        self.AddSimpleTool( self.CUSTOM_PAGE_ZOOM_ID, 'Zoom axes', loadIcon( "zoom-select.png" )
                           , ' Left mouse button to select rectangle zoom area\n Middle mouse button to pan axes\n Right click to zoom out.', wx.ITEM_CHECK )
        self.AddSimpleTool( self.CUSTOM_HOME_ID, 'Home', loadIcon( "zoom-original.png" )
                           , 'Restore the original view, i.e., reset zoom to fit image' )

        self.AddSeparator()
        self.AddSimpleTool( self.CUSTOM_SAVE_IMAGE_ID, 'Save image', loadIcon( "my-image-save.png" )
                            , 'Save plot contents to file' )
        #self.AddTool( self.CUSTOM_SAVE_IMAGE_ID, "Save image", loadIcon( "my-image-save.png" ), wx.NullBitmap
                        #, wx.ITEM_NORMAL, "Save image", "Save plot contents to file", None)

        #wx.EVT_TOOL( self, self.CUSTOM_PAN_LEFT_ID, self.onCustomPanLeft )
        #wx.EVT_TOOL( self, self.CUSTOM_PAN_RIGHT_ID, self.onCustomPanRight )
        wx.EVT_TOOL( self, self.CUSTOM_INC_CANVAS_WIDTH_ID, self.onIncCanvasWidth )
        wx.EVT_TOOL( self, self.CUSTOM_DEC_CANVAS_WIDTH_ID, self.onDecCanvasWidth )
        wx.EVT_TOOL( self, self.CUSTOM_FIT_CANVAS_ID, self.onFitCanvas )
        wx.EVT_TOOL( self, self.CUSTOM_ASPECT_AUTO_ID, self.onAspectAuto )
        wx.EVT_TOOL( self, self.CUSTOM_SAVE_IMAGE_ID, self.save )
        wx.EVT_TOOL( self, self.CUSTOM_PAGE_ZOOM_ID, self.onZoomButton )
        wx.EVT_TOOL( self, self.CUSTOM_HOME_ID, self.onHome )
        #self.Bind( wx.EVT_TEXT, self.onCanvasWidthCombo, id=self.CUSTOM_CANVAS_WIDTH_COMBO_ID )
        #self.canvasComboBox.Bind( wx.EVT_TEXT, self.onCanvasWidthCombo )

        self.mouseMoveID_ = None
        self.mousePressID_ = None
        self.mouseReleaseID_ = None
        self.buttonPressed_ = None
        self.rubberBander_ = MPLRubberBander( self.canvas )
        
        self.views_ = cbook.Stack()
        self.positions_ = cbook.Stack()

    def onZoomButton( self, event ):
        ''
        'Zoom button is pressed'
        ''
        if event.Checked():
            self.rubberBander_.start( zoomCallback = self.zoom, panCallback = self.pan )
        else:
            self.rubberBander_.stop( )

    def onHome( self, event ):
        ''
        ' Home button is pressed '
        ''
        lims = self.views_.home()
        if lims is not None:
            for i, a in enumerate( self.canvas.figure.get_axes() ):
                a.set_xlim( lims[ i ][ 0 ], lims[ i ][ 1 ] )
                a.set_ylim( lims[ i ][ 2 ], lims[ i ][ 3 ] )
                                
        self.panel.onZoomChanged( )
        self.panel.updateDrawOnIdle( )

    def zoom( self, selection ):
        ''
        ' zoom the current view'
        ''
        # selection = "coordPixel(left,right), coordAxes(left,right), axes
        
        # push the current view to define home if stack is empty
        if self.views_.empty(): self.push_current()
        
        leX = selection[ 1 ][ 0 ]
        riX = selection[ 1 ][ 2 ]
        leY = selection[ 1 ][ 1 ]
        riY = selection[ 1 ][ 3 ]
        selection[ 2 ].set_xlim( ( min( leX, riX ), max( leX, riX  ) ) )
        selection[ 2 ].set_ylim( ( min( leY, riY ), max( leY, riY  ) ) )
        
        # store new view sizes
        self.push_current()
        
        self.panel.onZoomChanged( )
        self.panel.updateDrawOnIdle( )

    def pan( self, pan ):
        ''
        ' pan the current view'
        ''
        # pan = "coordAxes(dx,dy), axes
        
        # push the current view to define home if stack is empty
        if self.views_.empty(): self.push_current()
        
        xlim = pan[ 1 ].get_xlim()
        ylim = pan[ 1 ].get_ylim()
        pan[ 1 ].set_xlim( ( xlim[ 0 ] - pan[0][0], xlim[ 1 ] - pan[0][0] ) );
        pan[ 1 ].set_ylim( ( ylim[ 0 ] - pan[0][1], ylim[ 1 ] - pan[0][1] ) );
        
        # store new view sizes
        self.push_current()

        self.panel.onPanChanged( )
        self.panel.updateDrawOnIdle( )

    def onIncCanvasWidth( self, event ):
#        print "onIncCanvasWidth", self.canvasComboBox.GetSelection( )
        curSelection = self.canvasComboBox.getSelection( )

        if curSelection < len( self.canvasWidthScales ) -1:
            if self.canvasWidthScales[ curSelection + 1 ] > self.canvasWidthScales[ curSelection ]:
                self.canvasComboBox.select( curSelection + 1 )
                self.onCanvasWidthCombo()
        elif curSelection == len( self.canvasWidthScales ) -1:
            self.canvasComboBox.select( 1 )
            self.onCanvasWidthCombo()

    def onDecCanvasWidth( self, event ):
#        print "onDecCanvasWidth", self.canvasComboBox.GetSelection( )
        curSelection = self.canvasComboBox.getSelection( )

        if curSelection > 0:
            if self.canvasWidthScales[ curSelection - 1 ] < self.canvasWidthScales[ curSelection ]:
                self.canvasComboBox.select( curSelection - 1 )
                self.onCanvasWidthCombo()

    def onCanvasWidthCombo( self, event = None ):
        curSelection = self.canvasComboBox.getSelection( )
 #       print "onCanvasWidthCombo: ", event, curSelection

        if curSelection == 0:
            self.EnableTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, False )
            self.EnableTool( self.CUSTOM_INC_CANVAS_WIDTH_ID, True )
        elif curSelection == len( self.canvasWidthScales ) -2:
            self.EnableTool( self.CUSTOM_INC_CANVAS_WIDTH_ID, False )
            self.EnableTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, True )
        elif curSelection == len( self.canvasWidthScales ) -1:
            self.EnableTool( self.CUSTOM_INC_CANVAS_WIDTH_ID, True )
            self.EnableTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, False )
        else:
            self.EnableTool( self.CUSTOM_INC_CANVAS_WIDTH_ID, True )
            self.EnableTool( self.CUSTOM_DEC_CANVAS_WIDTH_ID, True )

        newZoom = self.canvasWidthScales[ curSelection ];
        self.Refresh()

        if self.panel.canvasZoomWidth != newZoom:
            self.panel.canvasZoomWidth = newZoom
            self.panel.refresh()
            self.panel.SetupScrolling()

    def onFitCanvas( self, event ):
        self.canvasComboBox.SetValue( self.canvasWidthScalesNames[ -1 ] )
        self.onCanvasWidthCombo()

    def setAspectAuto( self, val = True ):
        ' Set auto aspect '
        if val:
            #if self.aspectTool.GetState() is not aui.AUI_BUTTON_STATE_CHECKED:
            self.aspectTool.SetState( aui.AUI_BUTTON_STATE_CHECKED )
            axis  = self.panel.figure.get_axes()
            for a in axis: a.set_aspect( 'auto' )
        else:
            #if self.aspectTool.GetState() is aui.AUI_BUTTON_STATE_CHECKED:
            self.aspectTool.SetState( aui.AUI_BUTTON_STATE_NORMAL )
            axis  = self.panel.figure.get_axes()
            for a in axis: a.set_aspect( 'equal' )

        self.Refresh()
        self.panel.needUpdateHack_ = True
        self.panel.updateDrawOnIdle()

    def onAspectAuto( self, event = None ):
        ' on auto aspect clicked '
        axis  = self.panel.figure.get_axes()
        for a in axis:
            if event.Checked():
                a.set_aspect( 'auto' )
            else:
                a.set_aspect( 'equal' )
        self.panel.needUpdateHack_ = True
        self.panel.updateDrawOnIdle()

    def save(self, evt):
        ' Fetch the required filename and file type. '
        filetypes, exts, filter_index = self.canvas._get_imagesave_wildcards()
        print("filter_index", filter_index, filetypes)
        default_file = "image"
        dlg = wx.FileDialog(self.panel, "Save to file", "", default_file,
                            filetypes,
                            wx.SAVE|wx.OVERWRITE_PROMPT|wx.CHANGE_DIR)
        dlg.SetFilterIndex(filter_index)
        if dlg.ShowModal() == wx.ID_OK:
            dirname  = dlg.GetDirectory()
            filename = dlg.GetFilename()

            format = exts[dlg.GetFilterIndex()]
            basename, ext = os.path.splitext(filename)
            print("ext = ", ext, format)
            if ext.startswith('.'):
                ext = ext[1:]
            if ext in ('svg', 'pdf', 'ps', 'eps', 'png') and format!=ext:
                #looks like they forgot to set the image type drop
                #down, going with the extension.
                warnings.warn('extension %s did not match the selected image type %s; going with %s'%(ext, format, ext), stacklevel=0)
                format = ext

            try:
                self.canvas.print_figure( os.path.join( dirname, basename + '.' + format), format=format)
            except Exception as e:
                error_msg_wx(str(e))
                
    def push_current(self):
        ''
        '    Push the current view limits and position onto the stack'
        ''
        lims = []; pos = []
        for a in self.canvas.figure.get_axes():
            xmin, xmax = a.get_xlim()
            ymin, ymax = a.get_ylim()
            lims.append( (xmin, xmax, ymin, ymax) )
            # Store both the original and modified positions
            pos.append( (
                    a.get_position(True).frozen(),
                    a.get_position().frozen() ) )
        self.views_.push(lims)
        self.positions_.push(pos)
        #self.set_history_buttons()
    # def push_current(...)

class wxMatplotPanelSimple( wx.Panel ):
    def __init__( self, renderPanel, color=None, dpi=None, **kwargs ):
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure
    # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, renderPanel, **kwargs )

        self.figure = Figure( None, dpi )
        #self.canvas = NoRepaintCanvas( self, -1, self.figure )

        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        sizer = wx.BoxSizer();
        sizer.Add( self.canvas, 1, wx.EXPAND )
        self.SetSizer( sizer )
        self.axes  = self.figure.add_subplot( 111 )
        self.axes.set_aspect( 'auto' )
        self.Bind(wx.EVT_SIZE, self._onSize)

    def _onSize( self, event = None ):
        pixels = tuple( [ self.GetSize()[0], self.GetSize()[1] ] )
        print(pixels)

        #if self.canvas.GetMinSize(  )[0] != pixels[0] or \
            #self.canvas.GetMinSize(  )[1] != pixels[1] :
        self.canvas.SetMinSize( pixels )

        self.figure.set_size_inches( float( pixels[ 0 ] )/self.figure.get_dpi(),
                                         float( pixels[ 1 ] )/self.figure.get_dpi() )
        self.canvas.draw()

class wxMatplotPanel( scrolled.ScrolledPanel  ):
    """
    The PlotPanel has a Figure and a Canvas.

    OnSize events simply set a flag, and the actually redrawing of the
    figure is triggered by an Idle event.
    """
    def __init__( self, renderPanel, color=None, dpi=None, **kwargs ):
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure

        # initialize Panel
        if 'id' not in list(kwargs.keys()):
            kwargs['id'] = wx.ID_ANY
        if 'style' not in list(kwargs.keys()):
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE

        scrolled.ScrolledPanel.__init__( self, renderPanel, **kwargs )
        self.renderPanel = renderPanel

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi )
        #self.canvas = NoRepaintCanvas( self, -1, self.figure )

        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )

        self.canvas.mpl_connect('button_press_event', self.onMousePress )
        self.canvas.mpl_connect('pick_event', self.onPick )

        sizer = wx.BoxSizer();
        sizer.Add( self.canvas, 1, wx.EXPAND )
        self.SetSizer( sizer )
        #self.SetAutoLayout(1)
        #self.SetupScrolling()

        self.SetColor( color )
        self._refresh = False
        self._updateDraw = False

        self.toolBar_ = None

        self.canvasZoomWidth = 1.0

        self.Bind(wx.EVT_IDLE, self._onIdle)
        self.Bind(wx.EVT_SIZE, self._onSize)

        self.resfreshCounter = 0
        self.needUpdateHack_ = False
        self.needDrawing = False
        self.refresh()

    def onPick( self, event ): pass

    def onMousePress( self, event ): pass
    
    def onZoomChanged( self ): pass

    def onPanChanged( self ): pass

    def getToolBar( self, parent = None ):
        if not self.toolBar_:
            self.toolBar_ = wxAUIMatplotPanelToolbar( self, self.canvas, parent )

        return self.toolBar_

    def SetColor( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()

        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )

    def _onSize( self, event ):
        self._refresh = True

    def _onIdle( self, evt ):
        if self.IsShownOnScreen():

            if self.needDrawing:
                self.redraw()

            if self._refresh:
                self.refresh()
                self._refresh = False
                
            if self._updateDraw:
                swatch = Stopwatch( True )

                self.canvas.draw()
                if self.needUpdateHack_:
                    self.needUpdateHack()
                    self.needUpdateHack_ = False
                self._updateDraw = False

                if self.canvasZoomWidth == 1.0:
                    self.SetupScrolling( False, False )
                print("draw: ",  swatch.duration())

    def updateDrawOnIdle( self ):
        self._updateDraw = True

    def resizeOnIdle( self ):
        self._refresh = True;

    def refresh( self ):
        swatch = Stopwatch( True )
        #pixels = tuple( self.GetParent().GetClientSize() )
        self.resfreshCounter += 1;

        pixels = tuple( [ int( self.GetSize()[0] * self.canvasZoomWidth ), int( self.GetSize()[1] * 1.0 ) ] )
        #    print self, self.resfreshCounter, pixels

        if self.canvas.GetMinSize(  )[0] != pixels[0] \
            or self.canvas.GetMinSize(  )[1] != pixels[1] :
            #print "resize canvas"
            #print "parent-size", self.renderPanel.GetSize()
            #print "self-size", self.GetSize()
            #print "tupel-size", pixels
        # to avoid _onSize loop under linux
        #if self.GetSize() != self.parent.GetClientSize():
        #if self.GetSize() != pixels:
            #self.SetSize( pixels )

        #self.canvas.SetSize( pixels )
            self.canvas.SetMinSize( pixels )

            self.figure.set_size_inches( float( pixels[ 0 ] )/self.figure.get_dpi(),
                                         float( pixels[ 1 ] )/self.figure.get_dpi() )

            adjust = True
            if hasattr( self, 'cbar' ):
                if self.cbar.active:
                    adjust = False;

            if pixels[ 0 ] > 50 and adjust:
                self.figure.subplotpars.update( left = 50.0 / pixels[ 0 ]
                                                , right = ( pixels[ 0 ] - 20.0 ) / pixels[ 0 ] )

                for a in self.figure.axes:
                    if hasattr( a, "update_params" ):
                        a.update_params()
                        #a.set_position( a.figbox, which = 'original' )
                        a.set_position( a.figbox, which = 'both' )

                #self.figure.subplots_adjust( left = 50.0 / pixels[ 0 ] )
                #self.figure.subplots_adjust( right = ( pixels[ 0 ] - 20.0 ) / pixels[ 0 ] )

        #self.canvas.draw()
        self.updateDrawOnIdle( )
        self.needUpdateHack_ = True;
           #print "refresh: ",  swatch.duration()

    def draw(self): pass # abstract, to be overridden by child classes

    def needUpdateHack( self ): pass
# END class wxMatplotPanel

class AppResourceWxMPL( AppResource, wxMatplotPanel ):
    """What is this??"""
    def __init__( self, parent, rendererSlot, propertyInspectorSlot ):
        AppResource.__init__( self, parent, rendererSlot, propertyInspectorSlot )
        wxMatplotPanel.__init__( self, rendererSlot, style = wx.NO_FULL_REPAINT_ON_RESIZE )

        self.SetColor( (255,255,255) );
        self.gci = None
        self.axes = None
        self.initAxes()

    def getRendererPanel( self ) : return self
    
    def createRendererPanel( self, parent): return self
    
    def initAxes( self ):
        #if not hasattr( self, 'axes' ):
            #self.axes = None

        #if self.axes is None:
        self.axes  = self.figure.add_subplot( 111 )
        self.axes.set_aspect( 'equal' )

    def redraw( self ):

        if isinstance( self.axes, list):
            for a in self.axes:
                self.figure.delaxes( a )
        else:
            self.figure.delaxes( self.axes )

        self.initAxes()
        self.draw()

    def draw( self ):
        ''
        ''
        ''
        if self.IsShownOnScreen():
            wx.BeginBusyCursor( wx.StockCursor( wx.CURSOR_WAIT ) )

            self.needDrawing = False
            
            for a in self.figure.axes:
                a.cla()
                
            self.drawData_()
            
            wx.EndBusyCursor( )
            self.updateDrawOnIdle()
        else:
            self.needDrawing = True

    def drawData_( self ):
        """Abstract interface, Define what we have to be drawn (needed from
        base class) is called while a draw event is fired."""
        Abstract_interface_need_to_be_derived