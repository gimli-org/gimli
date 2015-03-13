#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import wx
import wx.xrc as xrc

try:
    from agw import aui
    from agw.aui import aui_switcherdialog as ASD
except ImportError:  # if it's not there locally, try the wxPython lib.
    import wx.lib.agw.aui as aui
    from wx.lib.agw.aui import aui_switcherdialog as ASD

import wx.py as py
import wx.stc as stc
import wx.lib.scrolledpanel

import time

try:
    import pygimli as pg
    import pygimli.misc
except ImportError:
    sys.stderr.write(
        '''ERROR: cannot import the library 'pygimli'. Ensure that pygimli is in your PYTHONPATH \n''')
    sys.exit(1)
except Exception as e:
    print(e)

from pygimli.gui.resources import loadIcon
from pygimli.gui.controls import ResourceTree


class WorkSpace:

    def __init__(self):
        self.activeResource = None


def err(name):
    raise Exception(name)


class RedirectOutput:

    def __init__(self, frame, console, logFile=None):
        self.logFile = logFile
        self.frame = frame
        if console == "cout":  # sys.stdout
            self.style = 1  # black
        elif console == "cerr":  # sys.stderr
            self.style = 2  # red

        self.counter = 0

        self.timer = wx.PyTimer(self.writeOut)
        self.timer.Start(100)

        self.writeQueue_ = []
        self.frameInTab = False

    def flush(self):
        pass
        #import time
        # print "RedirectOutput flush", self.style, len(self.writeQueue_)
        # print self.writeQueue_

        #del self.writeQueue_[:]
        # while len(self.writeQueue_) > 0:
        # time.sleep(0.1)
        # pass

    def write(self, what):
        ##Do not put a print statement here!!##
        # self.logFile.write(what)
        self.writeQueue_.append(what)

    def writeOut(self, event=None):
        if len(self.writeQueue_):
            self.counter += 1

            if not self.frameInTab:
                tab = self.frame.GetParent()
                if tab:
                    tab.InsertPage(0, self.frame, self.frame.GetName())
                    tab.Layout()
                    # tab.Show()
                    self.frame.Show()
                    tab.SetSelection(0)
                    self.frameInTab = True

            if self.frame is not None:

                self.frame.SetReadOnly(False)
                self.frame.StartStyling(self.frame.GetCurrentPos(), 0xff)

                fi = None
                if self.logFile is not None:
                    fi = open(self.logFile, 'a+')

                while len(self.writeQueue_):
                    st = self.writeQueue_.pop(0)

                    #self.frame.AddText(str(("%d/%d-%d %s")%(counter,l,len(st), st)))
                    self.frame.AddText((str("%s") % st))
                    self.frame.EnsureCaretVisible()
                    self.frame.SetStyling(len(st), self.style)

                    if fi:
                        if self.style == 2:
                            fi.write('Error: ' + st.encode("utf-8"))
                        else:
                            fi.write(st.encode("utf-8"))

                self.frame.SetReadOnly(True)
                if fi:
                    fi.close()


class CustomStatusBar(wx.StatusBar):
    msgPosId_ = 0
    counterPosId_ = 1
    timerPosId_ = 2
    gauge_ = None
    statusMsg_ = None
    statusCounter_ = None

    def __init__(self, parent):
        wx.StatusBar.__init__(self, parent, -1)
        self.SetFieldsCount(3)
        self.SetStatusWidths([-1, -1, 200])

        self.gaugeFlag = False
        self.reposition()
        self.timer = wx.PyTimer(self.Notify)
        self.timer.Start(100)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_IDLE, self.OnIdle)

    def OnSize(self, event):
        # self.reposition()
        self.sizeChanged = True

    def OnIdle(self, event):
        if self.sizeChanged:
            self.reposition()

    def Notify(self):
        t = time.localtime(time.time())
        #st = time.strftime("%I:%M:%S", t)
        st = time.strftime("%d-%b-%Y   %I:%M:%S", t)
        self.SetStatusText(st, self.timerPosId_)

        # if self.statusMsg_:
        #st0 = ("%s" % self.statusMsg_)
        #self.SetStatusText(st0, self.msgPosId_)
        # else:
        #self.SetStatusText("Ready", self.msgPosId_)

        if self.statusCounter_:
            st0 = ("%d/%d" % self.statusCounter_())
            self.SetStatusText(st0, self.counterPosId_)
        else:
            self.SetStatusText(" ", self.counterPosId_)

        if self.gaugeFlag:
            self.gauge_.Show()
            self.gauge_.Pulse()
        else:
            self.gauge_.Hide()

    def pulsGauge(self, flag):
        self.gaugeFlag = flag

    def reposition(self):
        rect = self.GetFieldRect(1)
        if not self.gauge_:
            self.gauge_ = wx.Gauge(self, -1, 10, (0, 0), (10, rect.height))

        self.gauge_.SetSize((rect.width - 100, rect.height))
        self.gauge_.SetPosition((rect.x + 100, rect.y + 1))

        self.sizeChanged = False

    def setStatusMessage(self, msg, posID=0):
        #self.statusMsg_ = msg;
        self.SetStatusText(msg, posID)

    def setStatusCounter(self, counter):
        self.statusCounter_ = counter


class PyGIMLIMainFrame(wx.Frame):

    def __init__(self, ws):

        import locale
        if locale.localeconv()['decimal_point'] == ',':
            print("decimal point:", locale.localeconv()['decimal_point'])
            locale.setlocale(locale.LC_NUMERIC, 'C')

        self.ws = ws
        self.busyCursorWarning = False
        self.onIdleCmdQueue_ = list()

        self.xrc = xrc.EmptyXmlResource()
        xrcfile = 'pygimliMainFrame.xrc'

        if hasattr(sys, "frozen"):
            # frozen indicate windows.exe run??
            globPath = os.path.dirname(sys.argv[0])
        else:
            globPath = os.path.dirname(__file__)

        mainXRCFile = os.path.join(globPath, xrcfile)

        if not os.path.exists(mainXRCFile):
            raise IOError('Could not find xrc file "%s"; dying' % mainXRCFile)

        self.xrc.Load(mainXRCFile)

        self.PostCreate(self.xrc.LoadFrame(None, 'PyGIMLIMainFrame'))

        self.idleProzessPanel = self.xrc.LoadFrame(None, 'IdleProzessPanel')
        self.SetMinSize((800, 600))

        self.initFrameManager_()
        self.initDefaultMenuBar_()
        self.initToolBar_()
        self.initStatusBar_()

        self.initMainBINDS_()
        self.auiMgr.Update()

        self.stdout = None
        self.stderr = None

        self.plugins_ = dict()
        self.fileSuffixes = dict()

    def initFrameManager_(self):
        '''
            init the AUI-Manager and some default window targets
        '''
        self.auiMgr = aui.AuiManager()
        self.auiMgr.SetManagedWindow(self)

        self.auiMgr.AddPane(self.createRendererSlot_(self),
                            aui.AuiPaneInfo().Name("RenderWindow")
                            #                             .BestSize(wx.Size(250, 200)).MinSize(wx.Size(250, 200))
                            .CenterPane())

        self.auiMgr.AddPane(self.createResourceView_(self),
                            aui.AuiPaneInfo().Name(
                                "Resources").Caption("Resources")
                            .BestSize(wx.Size(300, 200)).MinSize(wx.Size(300, 200))
                            .CloseButton(True).MaximizeButton(True)
                            .Left().Layer(1).Position(0)
                            )

        self.auiMgr.AddPane(self.createPropertyInspectorSlot_(self),
                            aui.AuiPaneInfo().Name(
                                "Properties").Caption("Properties")
                            .BestSize(wx.Size(300, 200)).MinSize(wx.Size(300, 200))
                            .CloseButton(True).MaximizeButton(True)
                            .Left().Layer(1).Position(1)
                            )

        self.auiMgr.AddPane(self.createLogPane_(self),
                            aui.AuiPaneInfo().Name(
                                "Log and More").Caption("Log and More")
                            .BestSize(wx.Size(800, 150)).MinSize(wx.Size(500, 100))
                            .CloseButton(True).MaximizeButton(True)
                            .Bottom()
                            .Hide()
                            )

    def initDefaultMenuBar_(self):
        """Main menubar menues are come from mainXRCFile Init some default
        menubar entries File/Open, File/Quit."""

        # File/Open
        mbMenu = self.MenuBar.GetMenu(self.MenuBar.FindMenu("File"))

        item = wx.MenuItem(mbMenu, wx.NewId(), "&Open\tCtrl-O", "Open project")
        item.SetBitmap(loadIcon('document-open-16.png'))
        mbMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.onOpenFileDialog, id=item.GetId())

        # File/Quit
        item = wx.MenuItem(
            mbMenu,
            wx.NewId(),
            "&Quit\tCtrl-Q",
            "Quit application")
        item.SetBitmap(loadIcon('application-exit-16.png'))
        mbMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.onMainQuit, id=item.GetId())

        #item = wx.MenuItem(mbFileMenu, wx.NewId(), "Open &recent ...", "Open recent applications")
        # item.SetBitmap(loadIcon('application-exit-16.png'))
        # mbFileMenu.AppendItem(item)

    def initToolBar_(self):
        """Initilize some default toolbar entries."""
        self.toolBar_ = aui.AuiToolBar(self, -1, wx.DefaultPosition, wx.DefaultSize,
                                       aui.AUI_TB_DEFAULT_STYLE | aui.AUI_TB_PLAIN_BACKGROUND)
        #| aui.AUI_TB_OVERFLOW)

        self.toolBar_.SetToolBitmapSize(wx.Size(22, 22))
        tbquit = self.toolBar_.AddSimpleTool(wx.NewId(), "Quit", loadIcon("application-exit-16.png"), "Quit application"
                                             )
        wx.EVT_TOOL(self, tbquit.GetId(), self.onMainQuit)

        # self.toolBar_.Realize()
        self.auiMgr.AddPane(self.toolBar_, aui.AuiPaneInfo()
                            .Name("Main Toolbar").Caption("Main Toolbar")
                            .ToolbarPane().Top().Row(0).Position(0)
                            .LeftDockable(True).RightDockable(True))

    def initStatusBar_(self):
        """Initilize the default statusbar."""
        self.statusBar = CustomStatusBar(self)
        self.SetStatusBar(self.statusBar)

    def createRendererSlot_(self, parent):
        """"""
        self.rendererSlot = wx.Panel(parent)
        self.rendererSlot.SetSizer(wx.BoxSizer())
        return self.rendererSlot

    def createResourceView_(self, parent):
        """"""
        pane = ResourceTree(parent)
        self.resourceTree = pane
        return self.resourceTree

    def createPropertyInspectorSlot_(self, parent):
        """"""
        self.propertyInspectorSlot = wx.lib.scrolledpanel.ScrolledPanel(parent, -1, size=(14, 30),
                                                                        style=wx.TAB_TRAVERSAL | wx.SUNKEN_BORDER, name="PropertySlot")

        self.propertyInspectorSlot.SetSizer(wx.BoxSizer())
        self.propertyInspectorSlot.SetAutoLayout(1)
        self.propertyInspectorSlot.SetupScrolling()

        return self.propertyInspectorSlot

    def createLogPane_(self, parent):
        """"""
        pane = wx.Notebook(parent, -1, size=(810, 210), style=wx.BK_DEFAULT)

        self.log = stc.StyledTextCtrl(pane, -1)
        self.log.SetName("Log")
        self.log.SetUseHorizontalScrollBar(False)
        self.log.SetMarginWidth(1, 0)
        self.log.SetWrapMode(1)  # Turns on word wrap
        self.log.StyleClearAll()
        self.log.StyleSetSpec(1, "fore:BLACK")
        self.log.StyleSetSpec(2, "fore:RED")

        self.errLog = stc.StyledTextCtrl(pane, -1)
        self.errLog.SetName("Errors")
        self.errLog.SetUseHorizontalScrollBar(False)
        self.errLog.SetMarginWidth(1, 0)
        self.errLog.SetWrapMode(1)  # Turns on word wrap

        self.errLog.StyleClearAll()
        self.errLog.StyleSetSpec(1, "fore:BLACK")
        self.errLog.StyleSetSpec(2, "fore:RED")

        crust = wx.py.crust.Crust(pane, -1)
        self.log.Hide()
        self.errLog.Hide()

        #pane.AddPage(self.log, "Log")
        #pane.AddPage(self.errLog, "Errors")
        pane.AddPage(crust, "Crust")

        #pane.GetPage(0).Enable( False)

        self.logAndMore = pane
        return self.logAndMore

    def onClosePane(self, event):
        """"""
        window = event.GetPane().window

        if window == self.auiMgr.GetPane(self.logAndMore).window:
            self.GetMenuBar().Check(xrc.XRCID('mbViewLogAndMore'), False)
        elif window == self.auiMgr.GetPane(self.resourceTree).window:
            self.GetMenuBar().Check(xrc.XRCID('mbViewResources'), False)
        elif window == self.auiMgr.GetPane(self.propertyInspectorSlot).window:
            self.GetMenuBar().Check(xrc.XRCID('mbViewProperties'), False)

    def onSwitchViewPane(self, event):
        """"""
        if event.GetId() == xrc.XRCID('mbViewLogAndMore'):
            pane = self.auiMgr.GetPane(self.logAndMore)
        elif event.GetId() == xrc.XRCID('mbViewResources'):
            pane = self.auiMgr.GetPane(self.resourceTree)
        elif event.GetId() == xrc.XRCID('mbViewProperties'):
            pane = self.auiMgr.GetPane(self.propertyInspectorSlot)

        if event.IsChecked():
            pane.Show()
        else:
            pane.Hide()

        self.auiMgr.Update()

    def initMainBINDS_(self):
        '''
            Bind the default menu entities from XRC-file
        '''
        self.Bind(wx.EVT_CLOSE, self.onMainQuit)
        self.Bind(wx.EVT_MENU, self.onMainQuit, id=xrc.XRCID('tbQuit'))
        self.Bind(
            wx.EVT_MENU,
            self.onSwitchViewPane,
            id=xrc.XRCID('mbViewResources'))
        self.Bind(
            wx.EVT_MENU,
            self.onSwitchViewPane,
            id=xrc.XRCID('mbViewProperties'))
        self.Bind(
            wx.EVT_MENU,
            self.onSwitchViewPane,
            id=xrc.XRCID('mbViewLogAndMore'))
        self.Bind(wx.EVT_MENU, self.onAbout, id=xrc.XRCID('mbHelpAbout'))

        self.Bind(aui.EVT_AUI_PANE_CLOSE, self.onClosePane)
        self.Bind(wx.EVT_IDLE, self._onIdle)

    def redirectOutput(self, logFile=None):
        """Redirect stdout and stderr into a log console buffer."""
        sys.stdout = RedirectOutput(self.log, "cout", logFile)
        sys.stderr = RedirectOutput(self.errLog, "cerr", logFile)
        pass

    def addCommandToOnIdleQueue(self, cmd, args=[], label=""):
        if [cmd, args, label] not in self.onIdleCmdQueue_:
            self.onIdleCmdQueue_.append([cmd, args, label])

    def _onIdle(self, event):
        if len(self.onIdleCmdQueue_) > 0:
            # print len(self.onIdleCmdQueue_)
            gauge = xrc.XRCCTRL(self.idleProzessPanel, 'IdleGauge')

            if not self.idleProzessPanel.IsShown():
                gauge.SetRange(len(self.onIdleCmdQueue_) - 1)
                gauge.SetValue(0)
                self.idleProzessPanel.Layout()
                self.idleProzessPanel.Show()

                wx.BeginBusyCursor(wx.StockCursor(wx.CURSOR_WAIT))
            else:
                try:
                    gauge.SetValue(gauge.GetValue() + 1)
                except:
                    pass

            [cmd, args, name] = self.onIdleCmdQueue_[0]

            label = xrc.XRCCTRL(self.idleProzessPanel, 'IdleLabel')
            print(gauge.GetValue(), ": ", name)
            label.SetLabel("Prozessing: " + str(gauge.GetValue()) + "/" + str(gauge.GetRange())
                           + " ... " + name)

            try:
                if len(args) == 0:
                    print(name, cmd)
                    cmd()
                elif len(args) == 1:
                    cmd(args[0])
                elif len(args) == 2:
                    cmd(args[0], args[1])
            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print(e)

            self.onIdleCmdQueue_.pop(0)
        else:
            if self.idleProzessPanel.IsShown():
                wx.EndBusyCursor()
                self.idleProzessPanel.Hide()
            elif wx.IsBusy() and not self.busyCursorWarning:
                self.busyCursorWarning = True
                self.idleProzessPanel.Hide()
                wx.EndBusyCursor()
                err = wx.MessageDialog(
                    self,
                    'Hanging busy cursor found, probably something goes wrong. Please refer to the error log.',
                    'Something goes wrong.',
                    wx.OK | wx.ICON_WARNING)
                # err.ShowModal()
                if err.ShowModal() == wx.ID_OK:
                    self.busyCursorWarning = False

    def onOpenFileDialog(self, event=None):
        wildcard = str()
        for suffix in list(self.fileSuffixes.keys()):
            wildcard += self.fileSuffixes[suffix][0] + "|*" + suffix + "|"
        wildcard += "All files (*.*)|*.*"

        dlg = wx.FileDialog(self, message="Choose a file",
                            defaultDir=os.getcwd(),
                            defaultFile="",
                            wildcard=wildcard,
                            style=wx.OPEN | wx.CHANGE_DIR)  # | wx.FD_MULTIPLE)

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            print(('You selected %d files:' % len(paths)))

            self.openFile(paths[0])

    # Compare this with the debug above; did we change working dirs?
        print(("CWD: %s\n" % os.getcwd()))

        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()

    def openFile(self, path):
        (dirName, fileName) = os.path.split(path)
        (fileBaseName, fileExtension) = os.path.splitext(fileName)

        if fileExtension in self.fileSuffixes:
            print("Openfile: starting: ", self.fileSuffixes[fileExtension][1])
            app = self.createApplication(
                obj=self.fileSuffixes[fileExtension][1])
            self.fileSuffixes[fileExtension][2](app, path)
        else:
            err = wx.MessageDialog(self, 'File: ' + fileName + '\n'
                                   +
                                   'Cannot find an application associated to the file extension: '
                                   + fileExtension + '. Ignoring!', 'Something goes wrong while opening file.', wx.OK | wx.ICON_WARNING)
            err.ShowModal()

    def onMainQuit(self, event):
        # print "ok"
        #sys.stderr.write(" err\n")
        self.Destroy()
        sys.exit(0)

    def onAbout(self, event):
        """Show informations about pygimli."""
        from wx.lib.wordwrap import wordwrap
        info = wx.AboutDialogInfo()
        info.Name = "PyGI"
        info.Version = "0.9.0"
        info.Copyright = str(
            "(C) 2011 Carsten Rücker and Thomas Günther",
            'utf8')
        print(wx.PlatformInfo[1:])
        info.Description = wordwrap(
            "wyPython-" + wx.VERSION_STRING + " , ".join(wx.PlatformInfo[1:]) + ", \n" +
            " Running on python-" + sys.version.split()[0], 350, wx.ClientDC(self))
        info.WebSite = ("http://www.resistivity.net")
        info.Developers = [str("Carsten Rücker (carsten@resistivity.net)", 'utf8'),
                           str("Thomas Günther (thomas@resistivity.net)",
                               'utf8'),
                           ]

        info.License = wordwrap("licenseText", 500, wx.ClientDC(self))

        # Then we call wx.AboutBox giving it that info object
        wx.AboutBox(info)
        #-1-1 self.aboutGIMLiLabel.SetLabel(pg.version())
        # self.aboutGIMLiDialog.Show()

    def registerOpenFileSuffix(self, suffix, wildcard, cls, callback):
        """What is this?"""
        if suffix not in self.fileSuffixes:
            print("Register main open file suffix:", suffix, " for ", wildcard)
            self.fileSuffixes[suffix] = [wildcard, cls, callback]
        else:
            print(" there is already a definition for mainOpenFileSlot suffix: ", suffix,
                  "(" + self.fileSuffixes[suffix][0] + ")")

    def registerPlugins(self):
        """What is this?"""
        print("register plugins: ")
        pluginpath = os.path.dirname(__file__) + '/../apps/'
        paths = os.listdir(pluginpath)

        for p in paths:
            if not os.path.isdir(pluginpath + p) or p[0] is '.':
                # print p, "is not a plugin"
                continue

            pluginName = "pygimli.gui.apps." + p
            print("installing: ", pluginName)
            importCmd = "import " + pluginName + " as plugin"

            try:
                exec(importCmd)
            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print("import exception: ", e)
                continue

            if not hasattr(plugin, 'PluginApplication'):
                continue

            try:
                if hasattr(plugin, 'MainMenuBarNew_Item'):
                    menu = None

                    if self.MenuBar.FindMenu("&New") != -1:
                        menu = self.MenuBar.GetMenu(
                            self.MenuBar.FindMenu("&New"))
                    else:
                        menu = wx.Menu()
                        self.MenuBar.Insert(2, menu, "&New")

                    help = ""
                    if hasattr(plugin, 'MainMenuBarNew_ItemHelp'):
                        help = plugin.MainMenuBarNew_ItemHelp

                    item = wx.MenuItem(
                        menu,
                        wx.NewId(),
                        plugin.MainMenuBarNew_Item,
                        help)

                    # if bitmap is not None:
                    # item.SetBitmap(bitmap)

                    menu.AppendItem(item)
                    self.plugins_[item.GetId()] = plugin.PluginApplication
                    self.Bind(
                        wx.EVT_MENU,
                        self.createApplication,
                        id=item.GetId())

            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print("Exception in register MainMenuBar Entry", e)
                return

            try:
                if hasattr(plugin, 'MainOpenFileSuffix') and \
                        hasattr(plugin, 'MainOpenFileSlot') and \
                        hasattr(plugin, 'MainOpenWildcard'):

                    if isinstance(plugin.MainOpenFileSuffix, list):
                        for i, suffix in enumerate(plugin.MainOpenFileSuffix):
                            self.registerOpenFileSuffix(
                                suffix,
                                plugin.MainOpenWildcard[i],
                                plugin.PluginApplication,
                                plugin.MainOpenFileSlot)
                    else:
                        self.registerOpenFileSuffix(
                            plugin.MainOpenFileSuffix,
                            plugin.MainOpenWildcard,
                            plugin.PluginApplication,
                            plugin.MainOpenFileSlot)

            except Exception as e:
                import traceback
                traceback.print_exc(file=sys.stdout)
                print("Exception in register OpenFileSuffix", e)
                return

    def createApplication(self, event=None, obj=None):
        """What is this?"""

        if event is not None:
            obj = self.plugins_[event.GetId()]

        if obj is not None:
            app = obj(self, self.rendererSlot, self.propertyInspectorSlot)
            self.resourceTree.addItem(app, app.getName())
            self.resourceTree.selectItem(app)
            return app
        else:
            err("createApplication found no object to an application")

#from pygimli.utils import IPCServer, IPCThreadedTCPRequestHandler
#import threading

#import wx.lib.inspection
# wx.lib.inspection.InspectionTool().Show()


class PyGIMLIApp(wx.App):

    def __init__(self, options=None, args=None, ws=None):
        raise
        print("Mops")
        super().__init__(redirect=False)

        from optparse import OptionParser
        parser = OptionParser()
        parser.add_option(
            "",
            "--debug",
            dest="debug",
            action="store_true",
            help="Debug mode.",
            default=False)
        (options, args) = parser.parse_args()

        print(options, args)
        self.options = options
        self.args = args
        self.logFile = 'pygi.log'

        self.mainFrame = PyGIMLIMainFrame(ws=ws)

        if not options.debug:
            self.mainFrame.redirectOutput(self.logFile)

        self.SetTopWindow(self.mainFrame)
        self.mainFrame.Show()
        self.mainFrame.registerPlugins()

        #self.ipcServer_ = IPCServer(('localhost', 0), IPCThreadedTCPRequestHandler)
        #self.ipcServer_.environment = 'production'
        #ip, port = self.ipcServer_.server_address
        # print "starting ipc-server:", ip, port
        # Start a thread with the server -- that thread will then start one more thread for each request
        #server_thread = threading.Thread(target = self.ipcServer_.serve_forever)
        # Exit the server thread when the main thread terminates
        # server_thread.setDaemon(True)
        # server_thread.start()
        #ws.ipcServer = self.ipcServer_

    def start(self):
        for p in self.args:
            self.mainFrame.openFile(p)

        self.MainLoop()
