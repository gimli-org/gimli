# -*- coding: iso-8859-1 -*-

import wx

class ResourceTree( wx.TreeCtrl ):
    def __init__( self, parent ):
        self.parent = parent
        wx.TreeCtrl.__init__( self, parent )
        self.root = self.AddRoot( "Root" );
        self.SetPyData( self.root, parent )
        self.lastActiveItemData = None
        self.activeProfile = None

        self.Bind( wx.EVT_TREE_SEL_CHANGED, self.onSelectTreeObject )
        self.Bind( wx.EVT_CONTEXT_MENU, self.onPopupMenu )
        self.Bind( wx.EVT_KEY_DOWN, self.onKeyDown )
        
    def selectItem( self, data ):
        ''
        ''
        ''
        item  = None
        if hasattr( data, 'treeItem' ):
            item = data.treeItem
        else:
            item = data

        #print "selectItem", item
        if item:
            if item is not self.GetSelection():
                self.SelectItem( item )
        
    def addItem( self, data, name, parentNode = None ):
        if not parentNode:
            parentNode = self.root

        item = self.AppendItem( parentNode, name );
        self.SetPyData( item, data )
        data.treeItem = item
        self.Expand( parentNode )
        #self.UnselectAll( )
        
        return item
        
    def getSelectedData( self ):
        """
            Return the current selected item data
        """
        return self.GetPyData( self.GetSelection( ) )
        
    def onSelectTreeObject(self, event):
        try:
            itemData = event.GetEventObject().GetPyData( event.GetEventObject().GetSelection() )
        except:
            return
    
        if itemData is None or itemData is self.lastActiveItemData :
            return

        #print "onSelectTreeObject", itemData
        oldApplication = None
        newApplication = None
        
        if self.lastActiveItemData is not None:
            oldApplication = self.lastActiveItemData.parentResource
            if oldApplication is None:
                oldApplication = self.lastActiveItemData

            if hasattr( self.lastActiveItemData, 'activate' ):
                #print "self.lastActiveItemData.activate( False )", self.lastActiveItemData
                self.lastActiveItemData.activate( False )
                
        newApplication = itemData.parentResource
        if newApplication is None:
            newApplication = itemData
                
                
        # switch applications
        #print "old:", oldApplication
        #print "new:", newApplication
        if oldApplication != newApplication:
            if oldApplication is not None:
                oldApplication.activateApplication( False )
            if newApplication is not None:
                newApplication.activateApplication( True )
                
        self.lastActiveItemData = itemData
        #print "itemData.activate( True )", itemData
        itemData.activate( True )

        # store active resource in Workspace
        self.parent.ws.activeResource = itemData
        
        if itemData.parentResource is None:
            self.activeParent = itemData
        
        event.Skip()
        
    def onPopupMenu( self, evt ):
        
        menu = wx.Menu()
                
        self.popupDelete = wx.NewId()
        #self.popupReRead = wx.NewId()

        menu.AppendItem( wx.MenuItem(menu, self.popupDelete, "&Delete\tCTRL+D" ) )
        self.Bind(wx.EVT_MENU, self.onDeleteSelectedItem, id=self.popupDelete )
        
        # get local menus from application to add into popup menu
        app = self.GetPyData( self.GetSelection( ) )
        
        if app:
            if hasattr( app, 'menuSections' ):
                for i in range( len( app.menuSections ) ):
                    menuSection = app.menuSections[ i ]
                    menu.AppendSeparator(  )
                    for j in range( len( menuSection.items ) ):
                        name = menuSection.names[ j ]
                        func = menuSection.targetFunctions[ j ]
                        
                        item = wx.MenuItem( menu, wx.NewId(), name )
                        menu.AppendItem( item )
                        self.Bind( wx.EVT_MENU, func, id = item.GetId() )
            
        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        self.PopupMenu( menu )
        menu.Destroy()

    def onKeyDown(self, event ):
        #print "resourceTree::onKeyDown",event, event.GetKeyCode(), event.GetModifiers()
        
        if event.GetKeyCode() == wx.WXK_DELETE or \
            ( ( event.GetModifiers() == wx.MOD_CMD ) and event.GetKeyCode() == 68 ):

            self.onDeleteSelectedItem()
        
    def onDeleteSelectedItem( self, event = None ):
        item = self.GetPyData( self.GetSelection( ) )
        
        if not item:
            return
            
        try:
            # pls check the app will destroyed after deselect
            item.destroy( )
            prev = self.GetPrevSibling( self.GetSelection( ) )
            self.Delete( self.GetSelection( ) )
            self.SelectItem( prev )
        except:
            pass
