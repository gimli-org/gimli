#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os

import traceback

from PyQt4 import QtGui, QtCore

#import pygimli.gui.resources

# let ctrl-c on console abbort execution
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

from pygimli.gui.controls import ResourceTree

class WorkSpace:
    def __init__(self):
        self.activeResource = None

from spyderlib.widgets.internalshell import InternalShell
from code import InteractiveConsole
    
class PythonConsoleWidget(InternalShell):
    """
        Provides a custom widget to accept Python expressions.
    """
    def __init__(self, parent=None, mainApp=None):
        ns = {'win': self, 'app': mainApp, 'ws': mainApp.ws}
        super().__init__(parent, namespace=ns,
                         message="PyGIMLi internal console. \n\n"
                         "Access to application or workspace via: app, ws\n\n")
        #, message=("Spyder Internal Console\n\n" +
                                  #"This console is used to report application\n" +
                                  #"internal errors and to inspect Spyder\n" +
                                  #"internals with the following commands:\n" +
                                  #"  spy.app, spy.window, dir(spy)\n\n" +
                                  #"Please don't use it to run your code\n\n"))
        self.set_codecompletion_auto(True)
        self.set_calltips(True)
        self.interpreter.restore_stds()
        self.interpreter.push = self.push
    
    def push(self, line):
        """
            let all console output redirect to console .. everything else to the log
        """
        self.interpreter.redirect_stds()
        ret = InteractiveConsole.push(self.interpreter, "#coding=utf-8\n" + line)
        self.interpreter.restore_stds()
        return ret
    
class RedirectOutput():
    def __init__(self, frame, console, logFile=None):
        self.logFile = logFile
        self.frame = frame

        self._col = QtGui.QColor('black')
        self._style = 1

        if console == "cerr":
            self._col = QtGui.QColor('red')
            self._style = 2

        self.counter = 0

        self._timer = QtCore.QTimer()

        self._timer.connect(self._timer, QtCore.SIGNAL('timeout()'),
                            lambda: self.writeOut())
        self._timer.start(50)

        self.frame.insertPlainText("Start redirect for: " + console + "\n")
        self.frame.setReadOnly(True)
        self.writeQueue_ = []

    def write(self, what):
        ##Do not put a print statement here!!##
        self.writeQueue_.append(what)

    def writeOut(self):
        if len(self.writeQueue_):
            self.counter += 1

            if self.frame is not None:
                self.frame.setReadOnly(False)

                fi = None
                if self.logFile is not None:
                    fi = open(self.logFile, 'a+')

                while len(self.writeQueue_):
                    st = self.writeQueue_.pop(0)

                    self.frame.setTextColor(self._col)
                    #cursor = self.frame.textCursor()
                    #self.frame.textCursor().setPosition(QtGui.QTextCursor.end)
                    #cursor.setPosition(QtGui.QTextCursor.Down,
                                       #QtGui.QTextCursor.MoveAnchor)

                    #self.frame.setTextCursor(cursor)
                    #self.frame.textCursor().insertText(st)

                    self.frame.insertPlainText(st)
                    #self.frame.EnsureCaretVisible()

                    if fi:
                        if self._style == 2:
                            fi.write('Error: ' + st)
                        else:
                            fi.write(st)

                self.frame.setReadOnly(True)
                if fi:
                    fi.close()


class PyGUISystemMainFrame(QtGui.QMainWindow):
    """
    """
    def __init__(self, ws):
        """
        """
        super().__init__()

        import locale
        if locale.localeconv()['decimal_point'] == ',':
            print ("decimal point:", locale.localeconv()['decimal_point'])
            locale.setlocale(locale.LC_NUMERIC, 'C')

        self.ws = ws
        self._resourceTree = None

        self._createDefaultActions()

        self._initMenuBar()
        self._initToolBar()
        self._initStatusBar()
        self._createDockWindows()

        self._plugins = dict()
        self._fileSuffixes = dict()

        self._onIdleCmdQueue = []
        self._idleTimer = QtCore.QTimer()
        self.connect(self._idleTimer, QtCore.SIGNAL('timeout()'),
                     lambda: self._onIdle())
        self._idleTimer.start(100)

    def keyPressEvent(self, event):
        print("KEY", event.key())
        if type(event) == QtGui.QKeyEvent:
             #here accept the event and do something
            event.accept()
        else:
            event.ignore()

    def redirectOutput(self, logFile=None):
        """
            Redirect stdout and stderr into a log console buffer
        """
        sys.stdout = RedirectOutput(self._logAndMore, "cout", logFile)
        sys.stderr = RedirectOutput(self._logAndMore, "cerr", logFile)

    def registerOpenFileSuffix(self, suffix, wildcard, cls, callback):
        """
            What is this?
        """
        if suffix not in self._fileSuffixes:
            print("Register main open file suffix:", suffix, " for ", wildcard)
            self._fileSuffixes[suffix] = [wildcard, cls, callback]
        else:
            print("There is already a definition for fileSlot suffix: ",
                  suffix, "(" + self._fileSuffixes[suffix][0] + ")")

    def registerPlugins(self, pluginpath):
        """
            What is this?
        """
        import importlib
        paths = os.listdir(pluginpath)

        print("Registering plugins in: " + pluginpath)

        plugin = None
        newMB = None

        for p in paths:
            if not os.path.isdir(pluginpath + p) or p[0] is '.' or p[0] is '_':
                #print p, "is not a plugin"
                continue

            pluginName = "pygimli.gui.apps." + p
            print("Installing: ", pluginName)
            try:
                plugin = importlib.import_module(pluginName)
            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                print("Import pluging fails: ", e)
                continue

            if not hasattr(plugin, 'PluginApplication'):
                print("No valid PluginApplication found for:", pluginName)
                continue

            print("Setting menu entries for ", pluginName)
            try:
                if hasattr(plugin, 'MainMenuBarNew_Item'):

                    for mb in self.menuBar().findChildren(QtGui.QMenu):
                        if mb.title() == "&New":
                            newMB = mb
                    if newMB is None:
                        newMB = self.menuBar().addMenu('&New')

                    helpSt = ""
                    if hasattr(plugin, 'MainMenuBarNew_ItemHelp'):
                        helpSt = plugin.MainMenuBarNew_ItemHelp

                    action = newMB.addAction(plugin.MainMenuBarNew_Item)
                    action.setStatusTip(helpSt)
                    app = plugin.PluginApplication

                    # lambda in loops "scoping problem" in Python --
                    # the binding is late (lexical lookup at call-time)
                    # while you'd like it early (at def-time).
                    # resolve "fake default-value for argument"-idiom.
                    # lambda x=x: funct(x)as the

                    self.connect(action, QtCore.SIGNAL('triggered()'),
                                 lambda app=app: self.createApplication(app))

            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                print("Exception in register MainMenuBar Entry", e)
                continue

            try:
                if hasattr(plugin, 'MainOpenFileSuffix') and \
                   hasattr(plugin, 'MainOpenFileSlot') and \
                   hasattr(plugin, 'MainOpenWildcard'):

                    if type(plugin.MainOpenFileSuffix) is list:
                        for i, suffix in enumerate(plugin.MainOpenFileSuffix):
                            self.registerOpenFileSuffix(suffix,
                                plugin.MainOpenWildcard[i],
                                plugin.PluginApplication,
                                plugin.MainOpenFileSlot)
                    else:
                        self.registerOpenFileSuffix(plugin.MainOpenFileSuffix,
                                                    plugin.MainOpenWildcard,
                                                    plugin.PluginApplication,
                                                    plugin.MainOpenFileSlot)

            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                print("Exception in register OpenFileSuffix", e)
                return

            try:
                self.createApplication(plugin.PluginApplication)
            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                print(e)
                pass

    def createApplication(self, appDecl):
        """
            What is this?
        """
        if appDecl is not None:
            app = appDecl(self, self._rendererSlot, self._propertyView)
            self._resourceTree.addItem(app, app.name())
            self._resourceTree.selectItem(app)
            return app
        else:
            raise ("CreateApplication found no app to start on")

    def _createDefaultActions(self):
        """
        """
        self._exitAction = QtGui.QAction(
            QtGui.QIcon(':icons/application-exit'), 'Quit', self)
        self._exitAction.setShortcut('Ctrl+Q')
        self._exitAction.setStatusTip('Quit application')
        self.connect(self._exitAction,
                     QtCore.SIGNAL('triggered()'),
                     QtCore.SLOT('close()'))

        self._openAction = QtGui.QAction(
            QtGui.QIcon(':icons/application-open'), 'Open', self)
        self._openAction.setShortcut('Ctrl+O')
        self._openAction.setStatusTip('Open file')
        self.connect(self._openAction,
                     QtCore.SIGNAL('triggered()'), lambda: self._onMainOpen())

    def _initMenuBar(self):
        """
        """
        self.menuBar()
        self._fileMB = self.menuBar().addMenu('&File')
        self._viewMB = self.menuBar().addMenu('&View')
        self._fileMB.addAction(self._openAction)
        self._fileMB.addAction(self._exitAction)

    def _initToolBar(self):
        """
        """
        tb = self.addToolBar('File')
        tb.addAction(self._openAction)
        tb.addAction(self._exitAction)

    def _initStatusBar(self):
        """
        """
        self.statusBar()

    def _createDockWindows(self):
        """
        """
        self.setCentralWidget(self._createRendererSlot())

        resourcePane = QtGui.QDockWidget("Resources", self)
        resourcePane.setWidget(self._createResourceView())
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, resourcePane)
        self._viewMB.addAction(resourcePane.toggleViewAction())

        propertPane = QtGui.QDockWidget("Properties", self)
        propertPane.setWidget(self._createPropertyView())
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, propertPane)
        self._viewMB.addAction(propertPane.toggleViewAction())

        logPane = QtGui.QDockWidget("Logs and more", self)
        logPane.setWidget(self._createLogAndMore())
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, logPane)
        self._viewMB.addAction(logPane.toggleViewAction())

        consolePane = QtGui.QDockWidget("Console", self)
        consolePane.setWidget(PythonConsoleWidget(consolePane, self))
        self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, consolePane)
        self._viewMB.addAction(consolePane.toggleViewAction())

    def _createRendererSlot(self):
        """
        """
        #self._rendererSlot = QtGui.QMainWindow(self)
        #self._rendererSlot = QtGui.QMdiArea(self)
        self._rendererSlot = QtGui.QFrame(self)
        self._rendererSlot.setBaseSize(800, 600)
        #self._rendererSlot = QtGui.QMdiainWindow(self)
        #self._rendererSlot.setCentralWidget(QtGui.QWidget(self._rendererSlot))
        #self._rendererSlot.setSize(800,600)

        self._rendererSlot.setLayout(QtGui.QVBoxLayout())
        print("isDockNestingEnabled: ", self.isDockNestingEnabled())
        return self._rendererSlot

    def resourceTree(self):
        """
        """
        return self._createResourceView()

    def _createResourceView(self):
        """
        """
        if self._resourceTree is None:
            self._resourceTree = ResourceTree(self)
        return self._resourceTree

    def _createPropertyView(self):
        """
        """
        self._propertyView = QtGui.QListWidget()
        return self._propertyView

    def _createLogAndMore(self):
        """
        """
        self._logAndMore = QtGui.QTextEdit()
        return self._logAndMore

    def _onMainOpen(self):
        """
        """
        print("_onMainOpen(self):")

    def addCommandToOnIdleQueue(self, cmd, args=[], label=""):
        if [cmd, args, label] not in self._onIdleCmdQueue:
            self._onIdleCmdQueue.append([cmd, args, label])

    def _onIdle(self):
        """
        """
        if QtCore.QAbstractEventDispatcher.instance().hasPendingEvents() == True:
            return
        if len(self._onIdleCmdQueue) > 0:
            print("_onIdle")
            #gauge = xrc.XRCCTRL(self.idleProzessPanel, 'IdleGauge' )

            #if not self.idleProzessPanel.IsShown():
                #gauge.SetRange(len(self.onIdleCmdQueue_) -1)
                #gauge.SetValue(0)
                #self.idleProzessPanel.Layout()
                #self.idleProzessPanel.Show()

                #wx.BeginBusyCursor(wx.StockCursor(wx.CURSOR_WAIT))
            #else:
                #try:
                    #gauge.SetValue(gauge.GetValue() +1)
                #except:
                    #pass

            [cmd, args, name] = self._onIdleCmdQueue[0]

            #label = xrc.XRCCTRL(self.idleProzessPanel, 'IdleLabel' )
            #print(gauge.GetValue(), ": ", name)
            #label.SetLabel("Prozessing: " + str(gauge.GetValue()) +
                            #"/" + str(gauge.GetRange())+ " ... " + name)

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

            self._onIdleCmdQueue.pop(0)
        else:
            #if self.idleProzessPanel.IsShown():
                #wx.EndBusyCursor()
                #self.idleProzessPanel.Hide()
            #elif wx.IsBusy() and not self.busyCursorWarning:
                #self.busyCursorWarning = True
                #self.idleProzessPanel.Hide()
                #wx.EndBusyCursor()
                #err = wx.MessageDialog(self
                #, 'Hanging busy cursor found, probably something goes wrong.
                # Please refer to the error log.'
                #, 'Something goes wrong.'
                #, wx.OK | wx.ICON_WARNING)
                ##err.ShowModal()
                #if err.ShowModal() == wx.ID_OK:
                    #self.busyCursorWarning = False
            pass


class PyGIMLIApp(QtGui.QApplication):
    def __init__(self, options=None, args=None, ws=None):
        super().__init__(sys.argv)

        import argparse
        parser = argparse.ArgumentParser(description="GUI for Bert")
    
        parser.add_argument("--debug", dest="debug", action="store_true",
                            help="Debug mode.", default=False)
        args = parser.parse_args()

        self.args = args
        self.logFile = 'ogedit.log'

        self.mainFrame = PyGUISystemMainFrame(ws=ws)

        if not self.args.debug:
            self.mainFrame.redirectOutput(self.logFile)

        self.mainFrame.show()
        self.mainFrame.registerPlugins(os.path.dirname(__file__) +
                                       '/../apps/')
        self.mainFrame.resize(800, 600)

    def start(self):
        """
        """
        sys.exit(self.exec_())
