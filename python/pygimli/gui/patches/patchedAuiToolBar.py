# -*- coding: iso-8859-1 -*-

import wx

try:
    from agw import aui
    from agw.aui import aui_switcherdialog as ASD
    from agw.aui.aui_constants import *
except ImportError: # if it's not there locally, try the wxPython lib.
    import wx.lib.agw.aui as aui
    from wx.lib.agw.aui import aui_switcherdialog as ASD
    from wx.lib.agw.aui.aui_constants import *
    
class PatchedAuiToolBar( aui.AuiToolBar ):
    def __init__(self, parent, id=wx.ID_ANY, pos=wx.DefaultPosition,
                    size=wx.DefaultSize, style=AUI_TB_DEFAULT_STYLE):
        aui.AuiToolBar.__init__( self, parent, id, pos, size, style )
        
    def Realize(self):
        """
        Realizes the toolbar.

        This function should be called after you have added tools.
        """

        dc = wx.ClientDC(self)
        
        if not dc.IsOk():
            return False

        horizontal = True
        if self._style & AUI_TB_VERTICAL:
            horizontal = False

        # create the new sizer to add toolbar elements to
        sizer = wx.BoxSizer((horizontal and [wx.HORIZONTAL] or [wx.VERTICAL])[0])

        # add gripper area
        separator_size = self._art.GetElementSize(AUI_TBART_SEPARATOR_SIZE)
        gripper_size = self._art.GetElementSize(AUI_TBART_GRIPPER_SIZE)
        
        if gripper_size > 0 and self._gripper_visible:        
            if horizontal:
                self._gripper_sizer_item = sizer.Add((gripper_size, 1), 0, wx.EXPAND)
            else:
                self._gripper_sizer_item = sizer.Add((1, gripper_size), 0, wx.EXPAND)
        else:
            self._gripper_sizer_item = None
        
        # add "left" padding
        if self._left_padding > 0:
            if horizontal:
                sizer.Add((self._left_padding, 1))
            else:
                sizer.Add((1, self._left_padding))
        
        count = len(self._items)
        for i, item in enumerate(self._items):
        
            sizer_item = None
            kind = item.kind

            if kind == ITEM_LABEL:
                
                size = self._art.GetLabelSize(dc, self, item)
                sizer_item = sizer.Add((size.x + (self._tool_border_padding*2),
                                        size.y + (self._tool_border_padding*2)),
                                       item.proportion,
                                       item.alignment)
                if i+1 < count:
                    sizer.AddSpacer(self._tool_packing)
                

            elif kind in [ITEM_CHECK, ITEM_NORMAL, ITEM_RADIO]:
                
                size = self._art.GetToolSize(dc, self, item)
                sizer_item = sizer.Add((size.x + (self._tool_border_padding*2),
                                        size.y + (self._tool_border_padding*2)),
                                       0,
                                       item.alignment)
                # add tool packing
                if i+1 < count:
                    sizer.AddSpacer(self._tool_packing)

            elif kind == ITEM_SEPARATOR:
                
                if horizontal:
                    sizer_item = sizer.Add((separator_size, 1), 0, wx.EXPAND)
                else:
                    sizer_item = sizer.Add((1, separator_size), 0, wx.EXPAND)

                # add tool packing
                if i+1 < count:
                    sizer.AddSpacer(self._tool_packing)

            elif kind == ITEM_SPACER:
                
                if item.proportion > 0:
                    sizer_item = sizer.AddStretchSpacer(item.proportion)
                else:
                    sizer_item = sizer.Add((item.spacer_pixels, 1))
                    
            elif kind == ITEM_CONTROL:
                
                vert_sizer = wx.BoxSizer(wx.VERTICAL)
                vert_sizer.AddStretchSpacer(1)
                ctrl_sizer_item = vert_sizer.Add(item.window, 0, wx.EXPAND)
                vert_sizer.AddStretchSpacer(1)
                
                if self._style & AUI_TB_TEXT and \
                    self._tool_text_orientation == AUI_TBTOOL_TEXT_BOTTOM and \
                    item.GetLabel() != "":
                
                    s = self.GetLabelSize(item.GetLabel())
                    vert_sizer.Add((1, s.y))

                sizer_item = sizer.Add(vert_sizer, item.proportion, wx.EXPAND)
                min_size = item.min_size

                # proportional items will disappear from the toolbar if
                # their min width is not set to something really small
                if item.proportion != 0:
                    min_size.x = 1
                
               # if min_size.IsFullySpecified():
               #     sizer_item.SetMinSize(min_size)
               #     ctrl_sizer_item.SetMinSize(min_size)
                
                # add tool packing
                if i+1 < count:
                    sizer.AddSpacer(self._tool_packing)
                
            item.sizer_item = sizer_item
        

        # add "right" padding
        if self._right_padding > 0:
            if horizontal:
                sizer.Add((self._right_padding, 1))
            else:
                sizer.Add((1, self._right_padding))
        
        # add drop down area
        self._overflow_sizer_item = None

        if self._style & AUI_TB_OVERFLOW:
        
            overflow_size = self._art.GetElementSize(AUI_TBART_OVERFLOW_SIZE)
            if overflow_size > 0 and self._overflow_visible:
            
                if horizontal:
                    self._overflow_sizer_item = sizer.Add((overflow_size, 1), 0, wx.EXPAND)
                else:
                    self._overflow_sizer_item = sizer.Add((1, overflow_size), 0, wx.EXPAND)
            
            else:
            
                self._overflow_sizer_item = None
            
        # the outside sizer helps us apply the "top" and "bottom" padding
        outside_sizer = wx.BoxSizer((horizontal and [wx.VERTICAL] or [wx.HORIZONTAL])[0])

        # add "top" padding
        if self._top_padding > 0:
        
            if horizontal:
                outside_sizer.Add((1, self._top_padding))
            else:
                outside_sizer.Add((self._top_padding, 1))
        
        # add the sizer that contains all of the toolbar elements
        outside_sizer.Add(sizer, 1, wx.EXPAND)

        # add "bottom" padding
        if self._bottom_padding > 0:
        
            if horizontal:
                outside_sizer.Add((1, self._bottom_padding))
            else:
                outside_sizer.Add((self._bottom_padding, 1))

        del self._sizer # remove old sizer
        self._sizer = outside_sizer
        self.SetSizer(outside_sizer)

        # calculate the rock-bottom minimum size
        for item in self._items:
        
            if item.sizer_item and item.proportion > 0 and item.min_size.IsFullySpecified():
                item.sizer_item.SetMinSize((0, 0))
        
        self._absolute_min_size = self._sizer.GetMinSize()

        # reset the min sizes to what they were
        for item in self._items:
        
            if item.sizer_item and item.proportion > 0 and item.min_size.IsFullySpecified():
                item.sizer_item.SetMinSize(item.min_size)
        
        # set control size
        size = self._sizer.GetMinSize()
        self.SetMinSize(size)
        self._minWidth = size.x
        self._minHeight = size.y

        if self._style & AUI_TB_NO_AUTORESIZE == 0:
        
            cur_size = self.GetClientSize()
            new_size = self.GetMinSize()

            if new_size != cur_size:
            
                self.SetClientSize(new_size)
            
            else:
            
                self._sizer.SetDimension(0, 0, cur_size.x, cur_size.y)
            
        else:
        
            cur_size = self.GetClientSize()
            self._sizer.SetDimension(0, 0, cur_size.x, cur_size.y)
                    
        self.Refresh(False)
        return True
