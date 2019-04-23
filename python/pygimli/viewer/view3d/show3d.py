"""
Todo
----

+ log scale
"""


import matplotlib as mpl
import numpy as np
from shutil import copyfile
import sys
import time

import pygimli as pg

from PyQt5.QtCore import Qt, QPointF, QRect, QSize
from PyQt5.QtGui import (
    QPixmap, QPainter, QColor, QBrush, QIcon, QLinearGradient,
    QDoubleValidator
)
from PyQt5.QtWidgets import (
    QMainWindow, QFrame, QVBoxLayout, QToolBar, QComboBox, QPushButton,
    QFileDialog, QLineEdit, QWidget, QHBoxLayout, QSlider
)
vtki = pg.optImport('vtki', requiredFor="Proper visualization in 3D")


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

    def __init__(self, tmpMesh, parent=None):
        super(Show3D, self).__init__(parent)
        self.tmpMesh = tmpMesh
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

    def addMesh(self, mesh, cMap):
        """
        Add a mesh to the pyqt frame.

        Parameter
        ---------
        mesh: vtki.pointset.UnstructuredGrid
            The grid as it was read by vtki in vtkiview.
        cMap: str
            The MPL colormap that should be used to display parameters.
        """
        self.mesh = mesh
        self.vtk_widget.add_mesh(self.mesh, cmap=cMap)
        self.vtk_widget.add_bounds_axes()
        self.vtk_widget.reset_camera()

        # set the correctly chosen colormap
        if cMap.endswith('_r'):
            cMap = cMap[:-2]
            self.toolbar.btn_reverse.setChecked(2)
        # check if given cmap is in list and implement if it isn't
        if cMap not in CMAPS:
            self.toolbar.addExtraCMap(cMap)

        self.toolbar.cbbx_cmap.setCurrentText(cMap)
        self._allowSignals()

        _bounds = self.mesh.bounds
        self.toolbar.slice_x.setMinimum(_bounds[0])
        self.toolbar.slice_x.setMaximum(_bounds[1])
        self.toolbar.slice_x.setValue(0.5 * (_bounds[0] + _bounds[1]))
        self.toolbar.slice_y.setMinimum(_bounds[2])
        self.toolbar.slice_y.setMaximum(_bounds[3])
        self.toolbar.slice_y.setValue(0.5 * (_bounds[2] + _bounds[3]))
        self.toolbar.slice_z.setMinimum(_bounds[4])
        self.toolbar.slice_z.setMaximum(_bounds[5])
        self.toolbar.slice_z.setValue(0.5 * (_bounds[4] + _bounds[5]))

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
        # enable only when there is something to show
        self.toolbar.btn_apply.setEnabled(True)
        self.toolbar.btn_reset.setEnabled(True)
        self.toolbar.le_cmin.setEnabled(True)
        self.toolbar.le_cmax.setEnabled(True)

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

    def updateParameterView(self, param=None):
        """
        Change the view to given Parameter values.

        Parameter
        ---------
        param: Current text if the just triggered QComboBox

        Note
        ----
        Maybe overloaded.
        """
        if param is not None and param not in CMAPS and not isinstance(param, int):
            # change to the desired parameter distribution
            self.mesh.set_active_scalar(param)
            # update the minima and maxima
            self.toolbar.le_cmin.setText(self.extrema[param]['user']['min'])
            self.toolbar.le_cmax.setText(self.extrema[param]['user']['max'])

        cMap = self.toolbar.cbbx_cmap.currentText()
        if self.toolbar.btn_reverse.isChecked():
            cMap += '_r'

        if self.toolbar.btn_slice.isChecked():
            x_val = self.toolbar.slice_x.value()
            y_val = self.toolbar.slice_y.value()
            z_val = self.toolbar.slice_z.value()

            # multiBlock = self.mesh.slice_orthogonal(x=-0.05, y=-0.05, z=-0.05)
            # multiBlock.plot()
            # mesh = vtki.read(multiBlock)
            # self.vtk_widget.plot(multiBlock, cmap=cMap)
            # print(mesh)
        # else:
        mesh = self.mesh
        # print(dir(self.vtk_widget))
        self.vtk_widget.add_mesh(mesh, cmap=cMap)
        self.updateScalarBar()

    def updateScalarBar(self):
        """
        When user set limits are made and finished/accepted the color bar
        needs to change.
        """
        content_min = self.toolbar.le_cmin.text()
        content_max = self.toolbar.le_cmax.text()
        if content_min != '' or content_max != '':
            # get the user defined limits
            cmin = float(content_min)
            cmax = float(content_max)
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

    def exportMesh(self):
        f = QFileDialog.getSaveFileName(
            self, 'Export VTK', None, "VTK file (*.vtk)"
            )[0]
        if f:
            f = f + '.vtk' if not f.lower().endswith('.vtk') else f
            copyfile(self.tmpMesh, f)

    def resetExtrema(self):
        # get the active scalar/parameter that is displayed currently
        param = self.mesh.active_scalar_name
        self.extrema[param]['user']['min'] = self.extrema[param]['orig']['min']
        self.extrema[param]['user']['max'] = self.extrema[param]['orig']['max']

        # display correctly
        self.updateParameterView(param)

    def updateSlices(self):
        # slicer = self.sender()
        # orientation = slicer.toolTip()[4]
        # slider = getattr(self, 'slice_'.format(orientation.lower()))
        x_val = self.toolbar.slice_x.value()

    def enableSlicers(self):
        if self.toolbar.btn_slice.isChecked():
            self.toolbar.slice_x.setEnabled(True)
            self.toolbar.slice_y.setEnabled(True)
            self.toolbar.slice_z.setEnabled(True)
        else:
            self.toolbar.slice_x.setEnabled(False)
            self.toolbar.slice_x.setValue(self.toolbar.slice_x.minimum())
            self.toolbar.slice_y.setEnabled(False)
            self.toolbar.slice_y.setValue(self.toolbar.slice_y.minimum())
            self.toolbar.slice_z.setEnabled(False)
            self.toolbar.slice_z.setValue(self.toolbar.slice_z.minimum())

    def _allowSignals(self):
        # connect signals
        self.toolbar.cbbx_params.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.cbbx_cmap.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.btn_reverse.clicked.connect(self.updateParameterView)
        self.toolbar.btn_bbox.pressed.connect(self.toggleBbox)
        self.toolbar.btn_screenshot.clicked.connect(self.takeScreenShot)
        self.toolbar.btn_exportVTK.clicked.connect(self.exportMesh)
        self.toolbar.btn_apply.clicked.connect(self.updateScalarBar)
        self.toolbar.btn_reset.clicked.connect(self.resetExtrema)
        self.toolbar.btn_slice.clicked.connect(self.enableSlicers)
        self.toolbar.slice_x.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_y.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_z.sliderReleased.connect(self.updateParameterView)
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
        self.le_cmin.setEnabled(False)
        self.le_cmax = GLineEdit("The max of the current range")
        self.le_cmax.setEnabled(False)

        # if not hit by key his button accepts changes to the color range
        self.btn_apply = GButton(
            text="Apply",
            tooltip="Apply changes in color range"
            )
        self.btn_apply.setEnabled(False)

        # resets the color range to its original range
        self.btn_reset = GButton(
            text="Reset",
            tooltip="Reset changes in color range"
            )
        self.btn_reset.setEnabled(False)

        # slider for slicing
        self.slice_x = GSlider("The X location of the YZ slice")
        self.slice_x.setEnabled(False)
        self.slice_y = GSlider("The Y location of the XZ slice")
        self.slice_y.setEnabled(False)
        self.slice_z = GSlider("The Z location of the XY slice")
        self.slice_z.setEnabled(False)
        # and a botton to control them
        self.btn_slice = GButton(
            text="Slice",
            tooltip="Slice through volume",
            checkable=True
        )

        # button to take a screenshot
        self.btn_screenshot = GButton(
            text="Screenshot",
            tooltip="Save screenshot of the current scene"
            )

        self.btn_exportVTK = GButton(
            text="Export VTK",
            tooltip="Save displayed mesh as VTK"
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
        lt.addWidget(self.btn_slice)
        lt.addWidget(self.slice_x)
        lt.addWidget(self.slice_y)
        lt.addWidget(self.slice_z)
        lt.addStretch(1)
        lt.addWidget(self.btn_screenshot)
        lt.addWidget(self.btn_exportVTK)
        lt.setContentsMargins(0, 0, 0, 0)
        wt.setLayout(lt)
        self.addWidget(wt)

    def addExtraCMap(self, cMap):
        for icon, name in self._createPixmap([cMap]):
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


class _GDoubleSlider(QSlider):
    # https://gist.github.com/dennis-tra/994a65d6165a328d4eabaadbaedac2cc

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.decimals = 5
        self._max_int = 10 ** self.decimals

        super().setMinimum(0)
        super().setMaximum(self._max_int)

        self._min_value = 0.0
        self._max_value = 1.0

    @property
    def _value_range(self):
        return self._max_value - self._min_value

    def value(self):
        return float(super().value()) / self._max_int * self._value_range + self._min_value

    def setValue(self, value):
        super().setValue(int((value - self._min_value) / self._value_range * self._max_int))

    def setMinimum(self, value):
        if value > self._max_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._min_value = value
        self.setValue(self.value())

    def setMaximum(self, value):
        if value < self._min_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._max_value = value
        self.setValue(self.value())

    def minimum(self):
        return self._min_value

    def maximum(self):
        return self._max_value


class GSlider(_GDoubleSlider):

    def __init__(self, tooltip=None):
        super(GSlider, self).__init__(None)
        self.setOrientation(Qt.Horizontal)
        self.setToolTip(tooltip)
        # self.setMinimum(smin)
        # self.setMaximum(smax)
        # self.setSingleStep(1)


if __name__ == '__main__':
    app = Qt.QApplication(sys.argv)
    window = Show3D()
    sys.exit(app.exec_())
