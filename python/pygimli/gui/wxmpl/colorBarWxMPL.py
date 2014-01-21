# -*- coding: utf-8 -*-

import wx
import pygimli as g

from pygimli.gui.resources import loadXRC
from pygimli.gui.base import ManagedProperties

import matplotlib as mpl
from matplotlib.colors import Colormap

class ColorBarWxMPL( ManagedProperties ):
    def __init__( self, parent ):
        self.parent = parent
        ManagedProperties.__init__( self )
        self.cbar = None
        self.active = False
        
    def createPropertyPanel( self, panel = None ):
        xml = loadXRC( 'pygimli' )
        propertyPanel = None
        if panel:
            propertyPanel = xml.LoadPanel( panel, 'piColorBarWxMPL' )
            self.enablePanel = wx.xrc.XRCCTRL( panel, 'piColorBarPanel' )
            
            self.activeCheck = self.appendProperty( "cbActive"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarActiveCheckBox')
                                                    , ctrlEvent = wx.EVT_CHECKBOX
                                                    , targetFunct = self.activator )
            self.cbarTicks   = self.appendProperty( "cbNTicks"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarTicksSpinCtrl')
                                                    , ctrlEvent = wx.EVT_TEXT
                                                    , default = 5
                                                    , targetFunct = self.update )
            self.nColors     = self.appendProperty( "cbNColors"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarNColorsTextCtrl')
                                                    , ctrlEvent = wx.EVT_KILL_FOCUS
                                                    , default = 256
                                                    , valType = int
                                                    , targetFunct = self.update )
            self.nColors.setValidRange( 2, 256 )
            
            self.cmin       = self.appendProperty( "cbCmin"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarCminTextCtrl')
                                                    , ctrlEvent = wx.EVT_KILL_FOCUS
                                                    , valType = float
                                                    , targetFunct = self.update )
            self.cmax       = self.appendProperty( "cbCmax"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarCmaxTextCtrl')
                                                    , ctrlEvent = wx.EVT_KILL_FOCUS
                                                    , valType = float
                                                    , targetFunct = self.update )
            self.cAuto      = self.appendProperty( "cbAuto"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarAutoColorCheckBox')
                                                    , ctrlEvent = wx.EVT_CHECKBOX
                                                    , targetFunct = self.autoColor )
            self.logScale    = self.appendProperty( "cbLogScale"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarLogCheckBox')
                                                    , ctrlEvent = wx.EVT_CHECKBOX
                                                    , targetFunct = self.setLogScale )
            # better to set the following nativ through a propertygrid
            self.label       = self.appendProperty( "cbLabel"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarLabelTextCtrl')
                                                    , ctrlEvent = wx.EVT_KILL_FOCUS
                                                    , targetFunct = self.update )
            self.horizontal  = self.appendProperty( "cbHoriz" 
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarHorizRadioButton')
                                                    , ctrlEvent = wx.EVT_RADIOBUTTON
                                                    , targetFunct = self.update )
            self.vertical    = self.appendProperty( "cbVert"
                                                    , ctrl = wx.xrc.XRCCTRL( panel, 'piColorBarVertRadioButton')
                                                    , ctrlEvent = wx.EVT_RADIOBUTTON
                                                    , targetFunct = self.update )
        else:
            print("no panel given")

        return propertyPanel

    def setLogScale( self ):
        
        norm = None
        
        if self.logScale():
            if ( self.cmin() < 0 ):
                self.cmin.setVal( 1.0 )
            
            norm = mpl.colors.LogNorm()
        else:
            norm = mpl.colors.Normalize()
                
        if isinstance( self.parent.gci, list ):
            for gci in self.parent.gci:
                gci.set_norm( norm )
        else:
            self.parent.gci.set_norm( norm )

        self.update()
        
    def autoColor( self ):
        print("cBarAutoColor( self )", self.cAuto())
        if self.cbar:
            if self.cAuto():
                self.cmax.ctrl.Enable( False )
                self.cmin.ctrl.Enable( False )
                    
                cmin = 1e100
                cmax = -1e100
                if isinstance( self.parent.gci, list ):
                    for gci in self.parent.gci:
                        #gci.set_cmap( cmap )
                        #gci.set_clim( [ self.cmin(), self.cmax() ] )
                        #vals = ( 30, 50 )
                        cmin = min( cmin, gci.get_array().min() )
                        cmax = max( cmin, gci.get_array().max() )
                else:
                    vals = self.cbar.mappable.get_array( )
                    cmin = vals.min()
                    cmax = vals.max()

                self.cmin.ctrl.SetValue( "" )
                self.cmax.ctrl.SetValue( "" )
                self.cmin.val = cmin
                self.cmax.val = cmax
                
            else:
                self.cmax.ctrl.Enable( True )
                self.cmin.ctrl.Enable( True )
                self.cmin.ctrl.SetValue( str( self.cmin() ) )
                self.cmax.ctrl.SetValue( str( self.cmax() ) )

            self.update()
        
    def setGCI( self, gci ):
        if not self.cbar:
            self.cbar = g.mplviewer.createColorbar( gci )

    def activator( self, flag = None ):
        # possible bug: mehrfaches an und aus schalten verlangsamt alles.
       
        if flag is None:
            self.active = self.activeCheck()
        else:
            self.active = flag
            
        if self.active:
            try:
                self.enablePanel.Enable( True )
            except:
                pass
            
            if not self.cbar:
                if isinstance( self.parent.gci, list ):
                    gci = self.parent.gci[ -1 ]
                else:
                    gci = self.parent.gci
                    
                if gci.get_array( ) is not None:
                    if gci.get_array( ).min() != gci.get_array( ).max():
                        self.cbar = g.mplviewer.createColorbar( gci )
                        self.autoColor()
                        self.update()
                    else:
                        print("data range to small", gci.get_array( ).min(), gci.get_array( ).max())
                        self.activeCheck.setVal( False )
                        self.enablePanel.Enable( False )
                else:
                    print("no data asigned to mappable", self.parent.gci)
                    self.activeCheck.setVal( False )
                    self.enablePanel.Enable( False )
        else:
            if self.cbar:
                self.enablePanel.Enable( False )
                self.cbar.ax.cla()
                self.parent.figure.delaxes( self.cbar.ax )
                self.cbar = None;
                self.parent.redraw()
                
                #self.parent.updateDrawOnIdle()
                #self.parent.resizeOnIdle()

    def update( self ):
        if self.cbar:
            #print "horizontal", self.horizontal()
            nLevs = self.cbarTicks()
            self.cbar.set_label( self.label() )
            
            cmap = mpl.cm.get_cmap( 'jet',  self.nColors() )
            cmap.set_bad( [1.0, 1.0, 1.0, 0.0 ] )
            
            if isinstance( self.parent.gci, list ):
                for gci in self.parent.gci:
                    gci.set_cmap( cmap )
                    gci.set_clim( [ self.cmin(), self.cmax() ] )
            else:
                self.parent.gci.set_cmap( cmap )
                
            print(self.cmin(), self.cmax())
            if nLevs > 1:
                g.mplviewer.setCbarLevels( self.cbar, cMin = self.cmin(), cMax = self.cmax(), nLevs = self.cbarTicks() )

        self.parent.updateDrawOnIdle()
        
