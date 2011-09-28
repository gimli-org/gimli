#!/usr/bin/python
# -*- coding: utf-8 -*-

import wx

from listCtrlComboPopup import *
from checkListBoxComboPopup import *

class MainFrame( wx.Frame ):
    def __init__( self, *args, **kwds ):
        wx.Frame.__init__( self, *args, **kwds )

        self.sizer_  = wx.BoxSizer( wx.VERTICAL )

        self.listCtrlCombo = ListCtrlComboPopupControl( self
                                                , style=wx.CB_READONLY |wx.CB_DROPDOWN#|wx.NO_BORDER
                                                , size=(87,-1) )

        self.listCtrlCombo.addItems( ['item1', 'item2'] )
        self.sizer_.Add( self.listCtrlCombo, 0, wx.EXPAND )

        self.checkListBoxCombo = CheckListBoxComboPopupControl( self
                                                , style=wx.CB_READONLY |wx.CB_DROPDOWN#|wx.NO_BORDER
                                                , size=(87,-1) )

        self.checkListBoxCombo.addItems( ['item1', 'item2', 'item3', 'item4'] )
        self.sizer_.Add( self.checkListBoxCombo, 0, wx.EXPAND )

        self.SetSizer( self.sizer_ )


if __name__ == "__main__":
    testApp = wx.PySimpleApp( 0 )
    mainFrame = MainFrame( None, -1, "" )
    testApp.SetTopWindow(mainFrame)
    mainFrame.Show()
    testApp.MainLoop()