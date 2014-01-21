# -*- coding: utf-8 -*-
import wx
import wx.combo

class NullLog():
    def write( self, txt ):
        pass
        #print txt

class CheckListBoxComboPopupControl( wx.combo.ComboCtrl ):

    def __init__(self, *args, **kwargs):

        wx.combo.ComboCtrl.__init__( self, *args, **kwargs)
        self.popup = CheckListBoxComboPopup( )
        self.SetPopupControl( self.popup )

    def check( self, item ):
      #  print "CheckListBoxComboPopupControl.Select", item, self.popup.GetItemText( item )
        #self.SetValueWithEvent( self.popup.GetItemText( item ), True )
        pass

    def getChecked( self ):
        return self.popup.GetChecked()
        #self.SetValue( self.popup.GetItemText( item ) )

    def addItems( self, txt, checked = False ):
        for t in txt:
            self.popup.AddItem( t, checked )

    def SetValue( self, value ):
        """Give a tuple of strings for the ids of the checked boxes."""
        for i in self.popup.GetChecked():
            self.popup.Check( i, False )

        for v in value:
            if v == str:
                if v.isdigit():
                    self.popup.Check( int(v), True )
            else:
                self.popup.Check( v, True )

        wx.combo.ComboCtrl.SetValue( self, str( self.popup.GetChecked() ) )

    def GetValue( self ):
        """Return tuple of ids for the checked boxes."""
        wx.combo.ComboCtrl.SetValue( self, str( self.popup.GetChecked() ) )
        return self.popup.GetChecked()

class CheckListBoxComboPopup( wx.CheckListBox, wx.combo.ComboPopup):

    def __init__(self, log=None):
        self.dismissCallback = None
        if log:
            self.log = log
        else:
            self.log = NullLog()

        # Since we are using multiple inheritance, and don't know yet
        # which window is to be the parent, we'll do 2-phase create of
        # the ListCtrl instead, and call its Create method later in
        # our Create method.  (See Create below.)
        self.PostCreate( wx.PreCheckListBox() )

        # Also init the ComboPopup base class.
        wx.combo.ComboPopup.__init__(self)

    def AddItem(self, txt, checked = False ):
        self.Insert(txt, self.GetCount())
        self.Check( self.GetCount()-1, checked )

    def OnMotion(self, evt):
        pos = self.HitTest(evt.GetPosition())
        if pos >= 0:
            self.Select( pos )
            self.curitem = pos

    def OnLeftDown(self, evt):
        if not self.IsChecked( self.curitem ):
            self.Check( self.curitem, True )
        else:
            self.Check( self.curitem, False )

    # This is called immediately after construction finishes.  You can
    # use self.GetCombo if needed to get to the ComboCtrl instance.
    def Init(self):
        #self.log.write("CheckListBoxComboPopup.Init")
        self.value = -1
        self.curitem = -1

    # Create the popup child control.  Return true for success.
    def Create(self, parent):
        #self.log.write("CheckListBoxComboPopup.Create")
        wx.CheckListBox.Create( self, parent
                                #self, -1, (80, 50), wx.DefaultSize, sampleList
                            #              ,     style=wx.LC_LIST|wx.LC_SINGLE_SEL|wx.SIMPLE_BORDER
                                )

        self.Bind(wx.EVT_MOTION, self.OnMotion)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
        return True

    # Return the widget that is to be used for the popup
    def GetControl(self):
        #self.log.write("CheckListBoxComboPopup.GetControl")
        return self

    # Called immediately after the popup is shown
    def OnPopup(self):
        wx.combo.ComboPopup.OnPopup(self)

    # Called when popup is dismissed
    def OnDismiss(self):
        if self.dismissCallback:
            self.dismissCallback()
        wx.combo.ComboPopup.OnDismiss(self)

    # This is called to custom paint in the combo control itself
    # (ie. not the popup).  Default implementation draws value as
    # string.
    def PaintComboControl(self, dc, rect):
       # self.log.write("CheckListBoxComboPopup.PaintComboControl")
        wx.combo.ComboPopup.PaintComboControl(self, dc, rect)

    # Receives key events from the parent ComboCtrl.  Events not
    # handled should be skipped, as usual.
    def OnComboKeyEvent(self, event):
        wx.combo.ComboPopup.OnComboKeyEvent(self, event)

    # Implement if you need to support special action when user
    # double-clicks on the parent wxComboCtrl.
    def OnComboDoubleClick(self):
        wx.combo.ComboPopup.OnComboDoubleClick(self)

    # Return final size of popup. Called on every popup, just prior to OnPopup.
    # minWidth = preferred minimum width for window
    # prefHeight = preferred height. Only applies if > 0,
    # maxHeight = max height for window, as limited by screen size
    #   and should only be rounded down, if necessary.
    def GetAdjustedSize(self, minWidth, prefHeight, maxHeight):
        self.log.write("CheckListBoxComboPopup.GetAdjustedSize: %d, %d, %d" % (minWidth, prefHeight, maxHeight))
        try:
            prefHeight = self.GetCount() * ( self.GetItemHeight( ) + 7 )
        except Exception as e:
            print(e)            
            prefHeight = self.GetCount() * ( 17 + 7 )    
        return wx.combo.ComboPopup.GetAdjustedSize(self, minWidth, prefHeight, maxHeight)

    # Return true if you want delay the call to Create until the popup
    # is shown for the first time. It is more efficient, but note that
    # it is often more convenient to have the control created
    # immediately.
    # Default returns false.
    def LazyCreate(self):
        self.log.write("CheckListBoxComboPopup.LazyCreate")
        return wx.combo.ComboPopup.LazyCreate(self)
