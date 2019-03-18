"""
Todo
----

+ data-arg als pg.Vector und label-arg als name, dict optional
+ cmin/cmax
+ log scale
+ screenshot button
+ bounds on/off (able to siwtch on but not off)
"""


import matplotlib as mpl
import numpy as np
import sys
import time

import pygimli as pg

from PyQt5.QtCore import Qt, QPointF, QRect, QSize
from PyQt5.QtGui import (
    QPixmap, QPainter, QColor, QBrush, QIcon, QLinearGradient,
    QDoubleValidator
)
from PyQt5.QtWidgets import (
    QMainWindow, QFrame, QVBoxLayout, QToolBar, QComboBox,
    QPushButton, QFileDialog, QLineEdit, QWidget, QHBoxLayout
)
import vtki


# predefined color maps
CMAPS = ['viridis',
         'plasma',
         'Greys',
         'Blues',
         'Greens', 
         'Reds',
         'gray',
         'hot',
         'Spectral',
         'seismic',
         'Set3',
         'terrain',
         'rainbow',
         'jet'
         ]


class Show3D(QMainWindow):

    def __init__(self, parent=None):
        super(Show3D, self).__init__(parent)
        # storage for the minima and maxima
        self.extrema = {}
        self.setupWidget()

    def setupWidget(self):
        # create the frame
        self.frame = QFrame()
        vlayout = QVBoxLayout()
        vlayout.setContentsMargins(0, 0, 0, 0)

        # add the vtki interactor object
        self.vtk_widget = vtki.QtInteractor(self.frame)

        vlayout.addWidget(self.vtk_widget)

        self.frame.setLayout(vlayout)
        self.setupToolBar()

        self.setCentralWidget(self.frame)
        self.setWindowTitle("pyGIMLi 3D Viewer")

        self.show()

    def addMesh(self, mesh, cmap):
        """
        Add a mesh to the pyqt frame.

        Parameter
        ---------
        mesh: vtki.pointset.UnstructuredGrid
            The grid as it was read by vtki in vtkiview.
        cmap: str
            The MPL colormap that should be used to display parameters.
        """
        self.mesh = mesh
        self.vtk_widget.add_mesh(self.mesh, cmap=cmap)
        self.vtk_widget.add_bounds_axes()
        self.vtk_widget.reset_camera()

        # set the correctly chosen colormap
        if cmap.endswith('_r'):
            cmap = cmap[:-2]
            self.toolbar.btn_reverse.setChecked(2)
        # check if given cmap is in list and implement if it isn't
        if cmap not in CMAPS:
            self.toolbar.addExtraCMap(cmap)

        self.toolbar.cbbx_cmap.setCurrentText(cmap)
        self._allowSignals()

    def addDataToMesh(self, data):
        """
        Supply data to visualize.

        Parameter
        ---------
        data: dict
            Dictionary of cell values, sorted by key.

        Note
        ----
        This apparantly needs to happen when using the gui since on call the
        cell_arrays will be emptied...
        """
        for k, v in data.items():
            self.mesh._add_cell_scalar(v, k)
            self.extrema[k] = {
                'orig': {'min': str(min(v)), 'max': str(max(v))},
                'user': {'min': str(min(v)), 'max': str(max(v))}
                }
        # supply the combobox with the names to choose from for display
        self.toolbar.cbbx_params.addItems(self.mesh.scalar_names)
        # set the first cMin/cMax
        self.toolbar.le_cmin.setText(self.extrema['_Attribute']['orig']['min'])
        self.toolbar.le_cmax.setText(self.extrema['_Attribute']['orig']['max'])

    def setupToolBar(self):
        """
        Set up the toolbar and provide actions.
        """
        # init toolbar
        self.toolbar = GToolBar()
        # add to mainwindow
        self.addToolBar(self.toolbar)

    def updateParameterView(self, param):
        """
        Change the view to given Parameter values.

        Parameter
        ---------
        param: Current text if the just triggered QComboBox

        Note
        ----
        Maybe overloaded.
        """
        if param not in CMAPS and not isinstance(param, int):
            # change to the desired parameter distribution
            self.mesh.set_active_scalar(param)
            # update the minima and maxima
            self.toolbar.le_cmin.setText(self.extrema[param]['user']['min'])
            self.toolbar.le_cmax.setText(self.extrema[param]['user']['max'])

        cmap = self.toolbar.cbbx_cmap.currentText()
        if self.toolbar.btn_reverse.isChecked():
            cmap += '_r'
        self.vtk_widget.add_mesh(self.mesh, cmap=cmap)
        # self.vtk_widget.update()
        self.updateScalarBar()

    def updateScalarBar(self):
        """
        When user set limits are made and finished/accepted the color bar
        needs to change.
        """
        # get the user defined limits
        cmin = float(self.toolbar.le_cmin.text())
        cmax = float(self.toolbar.le_cmax.text())
        if cmax > cmin:
            # get the active scalar/parameter that is displayed currently
            param = self.mesh.active_scalar_name

            # update the user extrema
            self.extrema[param]['user']['min'] = str(cmin)
            self.extrema[param]['user']['max'] = str(cmax)

            # actually update the range
            self.vtk_widget.update_scalar_bar_range([cmin, cmax], name=param)
        self.vtk_widget.update()

    def _checkDecimalPoint(self):
        self.toolbar.le_cmin.setText(self.toolbar.le_cmin.text().replace(',', '.'))
        self.toolbar.le_cmax.setText(self.toolbar.le_cmax.text().replace(',', '.'))

    def toggleBbox(self):
        checked = not self.toolbar.btn_bbox.isChecked()
        self.vtk_widget.add_bounds_axes(
            show_xaxis=checked,
            show_yaxis=checked,
            show_zaxis=checked,
            show_xlabels=checked,
            show_ylabels=checked,
            show_zlabels=checked
        )
        self.vtk_widget.update()

    def takeScreenShot(self):
        fname = QFileDialog.getSaveFileName(
            self, 'Open File', None, "Image files (*.jpg *.png)"
            )[0]
        if fname:
            if not len(fname.split('.')) == 2:
                fname += '.png'
            self.vtk_widget.screenshot(fname)

    def resetExtrema(self):
        # get the active scalar/parameter that is displayed currently
        param = self.mesh.active_scalar_name
        self.extrema[param]['user']['min'] = self.extrema[param]['orig']['min']
        self.extrema[param]['user']['max'] = self.extrema[param]['orig']['max']

        # display correctly
        self.updateParameterView(param)

    def _allowSignals(self):
        # connect signals
        self.toolbar.cbbx_params.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.cbbx_cmap.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.btn_reverse.clicked.connect(self.updateParameterView)
        self.toolbar.btn_bbox.pressed.connect(self.toggleBbox)
        self.toolbar.btn_screenshot.clicked.connect(self.takeScreenShot)
        self.toolbar.btn_apply.clicked.connect(self.updateScalarBar)
        self.toolbar.btn_reset.clicked.connect(self.resetExtrema)
        self.toolbar.le_cmin.editingFinished.connect(self.updateScalarBar)
        self.toolbar.le_cmax.editingFinished.connect(self.updateScalarBar)


class GToolBar(QToolBar):
    """
    Provide the toolbar for the 3D viewer with all buttons.
    This just provides the graphical features.
    """

    def __init__(self, parent=None):
        super(GToolBar, self).__init__(parent)
        self.setupWidget()

    def setupWidget(self):
        """
        Get some actions..
        """
        # restrict from shifting the toolbar since the wider buttons would
        # make an awful wide sidebar
        self.setMovable(False)
        # combobox to choose the given parameter from
        self.cbbx_params = GComboBox("Scroll/Click to select parameter")

        # combobox to choose the colormap
        self.cbbx_cmap = GComboBox("Scroll/Click to select color map")
        for icon, name in self._createPixmap():
            self.cbbx_cmap.addItem(icon, name)
            self.cbbx_cmap.setIconSize(QSize(40, 15))

        # checkbox to reverse the chosen color scheme
        self.btn_reverse = GButton(
            text="_r",
            tooltip="Reverse the chosen color map",
            checkable=True,
            size=[24, 24]
        )

        # button to make bounding box visible
        self.btn_bbox = GButton(
            text="BBox",
            tooltip="Toggle data axis grid",
            checkable=True
        )
        self.btn_bbox.setChecked(True)

        # cMin and cMax
        self.le_cmin = GLineEdit("The min of the current range")
        self.le_cmax = GLineEdit("The max of the current range")

        # if not hit by key his button accepts changes to the color range
        self.btn_apply = GButton(
            text="Apply",
            tooltip="Apply changes in color range"
            )
        # resets the color range to its original range
        self.btn_reset = GButton(
            text="Reset",
            tooltip="Reset changes in color range"
            )

        # button to take a screenshot
        self.btn_screenshot = GButton(
            text="screenshot",
            tooltip="Save screenshot of the scene"
            )

        # widget for the toolbar to better organize the widgets
        wt = QWidget()
        lt = QHBoxLayout()
        lt.addWidget(self.cbbx_params)
        lt.addWidget(self.cbbx_cmap)
        lt.addWidget(self.btn_reverse)
        lt.addWidget(self.btn_bbox)
        lt.addWidget(self.le_cmin)
        lt.addWidget(self.le_cmax)
        lt.addWidget(self.btn_apply)
        lt.addWidget(self.btn_reset)
        lt.addStretch(1)
        lt.addWidget(self.btn_screenshot)
        lt.setContentsMargins(0, 0, 0, 0)
        wt.setLayout(lt)
        self.addWidget(wt)

    def addExtraCMap(self, cmap):
        for icon, name in self._createPixmap([cmap]):
            self.cbbx_cmap.addItem(icon, name)
            self.cbbx_cmap.setIconSize(QSize(40, 15))

    def _createPixmap(self, cmaps=CMAPS):
        """
        Generator to create icons shown in the comboboxes where the colormap
        is selected.
        """
        w = 60
        h = 20
        for cmap in cmaps:
            cMap = mpl.cm.get_cmap(cmap)

            px = QPixmap(w, h)
            p = QPainter(px)
            gradient = QLinearGradient(QPointF(0, 0), QPointF(w, 0))

            for t in np.linspace(0., 1., 15):
                rgba = cMap(t)
                gradient.setColorAt(t, QColor(rgba[0]*255, rgba[1]*255, rgba[2]*255))

            brush = QBrush(gradient)
            p.fillRect(QRect(0, 0, w, h), brush)
            p.end()

            yield QIcon(px), cmap


class GButton(QPushButton):

    def __init__(self, text=None, tooltip=None, checkable=False, size=None):
        super(GButton, self).__init__(None)
        self.setText(text)
        self.setToolTip(tooltip)
        self.setCheckable(checkable)
        if size is not None:
            self.setFixedSize(*size)


class GLineEdit(QLineEdit):

    def __init__(self, tooltip=None):
        super(GLineEdit, self).__init__(None)
        self.setToolTip(tooltip)
        self.setFixedWidth(150)
        # restrict acceptance to numbers only in that range
        # NOTE: will accept numbers outside of that range but the the
        # 'editingFinished' signal (hit enter/return when finished) won't
        self.setValidator(QDoubleValidator(-1e9, 1e9, 2, self))


class GComboBox(QComboBox):

    def __init__(self, tooltip=None):
        super(GComboBox, self).__init__(None)
        self.setToolTip(tooltip)


if __name__ == '__main__':
    app = Qt.QApplication(sys.argv)
    window = Show3D()
    sys.exit(app.exec_())
