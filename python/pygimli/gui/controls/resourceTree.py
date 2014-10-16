# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore

class ResourceTreeItem(QtGui.QTreeWidgetItem):
    def __init__(self, parent, name):
        QtGui.QTreeWidget.__init__(self, parent, name)
        self._data = None
        
    def setPyData(self, data):
        self.data_ = data
        
    def pyData(self):
        return self.data_

class ResourceTree(QtGui.QTreeWidget):
    
    def __init__(self, parent):
        self.parent = parent
        QtGui.QTreeWidget.__init__(self, parent)
        
        #self.setTitle('Resources')
        self.appData_ = dict()
        #self.root = self.addItem(None, "Root")
        
        #self.SetPyData(self.root, parent)
        self.lastActiveItemData = None
        self.activeProfile = None

        
        self.itemClicked.connect(lambda: self.onSelectTreeObject())
            
        #self.Bind(wx.EVT_TREE_SEL_CHANGED, self.onSelectTreeObject)
        #self.Bind(wx.EVT_CONTEXT_MENU, self.onPopupMenu)
        #self.Bind(wx.EVT_KEY_DOWN, self.onKeyDown)
    
    def getSelectedData(self):
        """
            Return the current selected item data.
        """
        return self.currentItem().pyData()
    
    def selectItem(self, data):
        """
        """
        item  = None
        if hasattr(data, 'treeItem'):
            item = data.treeItem
        else:
            item = data

        if item:
            if item is not self.currentItem():
                self.setCurrentItem(item)
        
    def addItem(self, data, name, parentNode=None):
        """
        """
        item = None
        
        if not parentNode:
            parentNode = self.invisibleRootItem()

        item = ResourceTreeItem(parentNode, name)
        item.setText(0, name)
        parentNode.addChild(item)
        
        if data is not None:
            item.setPyData(data)
            data.treeItem = item
        
        self.expandItem(parentNode)
        #self.UnselectAll()
        
        return item
        
    def onSelectTreeObject(self, item=None):
        """
            
        """
        print("onSelectTreeObject(self, item=None):")
        
        try:
            itemData = self.getSelectedData()
        except:
            return
    
        print("Data:", itemData) 
    
        if itemData is None or itemData is self.lastActiveItemData :
            return

        #print "onSelectTreeObject", itemData
        oldApplication = None
        newApplication = None
        
        if self.lastActiveItemData is not None:
            oldApplication = self.lastActiveItemData.parentResource
            if oldApplication is None:
                oldApplication = self.lastActiveItemData

            if hasattr(self.lastActiveItemData, 'activate'):
                #print "self.lastActiveItemData.activate(False)", self.lastActiveItemData
                self.lastActiveItemData.activate(False)
                
        newApplication = itemData.parentResource
        if newApplication is None:
            newApplication = itemData
                
                
        # switch applications
        #print "old:", oldApplication
        #print "new:", newApplication
        if oldApplication != newApplication:
            if oldApplication is not None:
                oldApplication.activateApplication(False)
            if newApplication is not None:
                newApplication.activateApplication(True)
                
        self.lastActiveItemData = itemData
        #print "itemData.activate(True)", itemData
        itemData.activate(True)

        # store active resource in Workspace
        self.parent.ws.activeResource = itemData
        
        if itemData.parentResource is None:
            self.activeParent = itemData
        
    def onPopupMenu(self, evt):
        
        menu = wx.Menu()
                
        self.popupDelete = wx.NewId()
        #self.popupReRead = wx.NewId()

        menu.AppendItem(wx.MenuItem(menu, self.popupDelete, "&Delete\tCTRL+D"))
        self.Bind(wx.EVT_MENU, self.onDeleteSelectedItem, id=self.popupDelete)
        
        # get local menus from application to add into popup menu
        app = self.GetPyData(self.GetSelection())
        
        if app:
            if hasattr(app, 'menuSections'):
                for i in range(len(app.menuSections)):
                    menuSection = app.menuSections[ i ]
                    menu.AppendSeparator()
                    for j in range(len(menuSection.items)):
                        name = menuSection.names[ j ]
                        func = menuSection.targetFunctions[ j ]
                        
                        item = wx.MenuItem(menu, wx.NewId(), name)
                        menu.AppendItem(item)
                        self.Bind(wx.EVT_MENU, func, id = item.GetId())
            
        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        self.PopupMenu(menu)
        menu.Destroy()

    def onKeyDown(self, event):
        #print "resourceTree::onKeyDown",event, event.GetKeyCode(), event.GetModifiers()
        
        if event.GetKeyCode() == wx.WXK_DELETE or \
            ((event.GetModifiers() == wx.MOD_CMD) and event.GetKeyCode() == 68):

            self.onDeleteSelectedItem()
        
    def onDeleteSelectedItem(self, event = None):
        item = self.GetPyData(self.GetSelection())
        
        if not item:
            return
            
        try:
            # pls check the app will destroyed after deselect
            item.destroy()
            prev = self.GetPrevSibling(self.GetSelection())
            self.Delete(self.GetSelection())
            self.SelectItem(prev)
        except:
            pass
