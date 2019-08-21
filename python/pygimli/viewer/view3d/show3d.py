"""
Todo
----

+ cache settings for next run
+ log scale
+ coverage, i.e. with opacity?
+ slider to look in volume
"""


import matplotlib as mpl
import numpy as np
from shutil import copyfile
import sys

import pygimli as pg

from PyQt5.QtCore import Qt, QPointF, QRect, QSize
from PyQt5.QtGui import (
    QPixmap, QPainter, QColor, QBrush, QIcon, QLinearGradient,
    QDoubleValidator
)
from PyQt5.QtWidgets import (
    QMainWindow, QFrame, QVBoxLayout, QToolBar, QComboBox, QPushButton,
    QFileDialog, QLineEdit, QWidget, QHBoxLayout, QSlider, QSplitter,
    QGroupBox, QLabel, QLineEdit
)
pyvista = pg.optImport('pyvista', requiredFor="proper visualization in 3D")


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

        # add the pyvista interactor object
        self.pyvista_widget = pyvista.QtInteractor(self.frame)

        vlayout = QVBoxLayout()
        vlayout.setContentsMargins(0, 0, 0, 0)
        # vlayout.addWidget(self.pyvista_widget)

        self.toolbar = GToolBar()

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.toolbar)
        splitter.addWidget(self.pyvista_widget)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 5)

        vlayout.addWidget(splitter)

        self.frame.setLayout(vlayout)
        # self.setupToolBar()

        self.setCentralWidget(self.frame)
        self.setWindowTitle("pyGIMLi 3D Viewer")

        self.show()

    def addMesh(self, mesh, cMap):
        """
        Add a mesh to the pyqt frame.

        Parameter
        ---------
        mesh: pyvista.pointset.UnstructuredGrid
            The grid as it was read by pyvista in vistaview.
        cMap: str
            The MPL colormap that should be used to display parameters.
        """
        self.mesh = mesh
        self._actor = self.pyvista_widget.add_mesh(self.mesh, cmap=cMap)
        self.pyvista_widget.show_bounds()
        self.pyvista_widget.reset_camera()

        # set the correctly chosen colormap
        if cMap.endswith('_r'):
            cMap = cMap[:-2]
            self.toolbar.btn_reverse.setChecked(2)
        # check if given cmap is in list and implement if it isn't
        if cMap not in CMAPS:
            self.toolbar.addExtraCMap(cMap)

        self.toolbar.cbbx_cmap.setCurrentText(cMap)
        self.allowMeshParameters()
        self._allowSignals()

        # set slicers to center after they're enabled
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

    def allowMeshParameters(self):
        """
        Make data from the given mesh accessible via GUI.

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

        # FIXME: what about the point arrays?!
        for k, v in self.mesh.cell_arrays.items():
            self.mesh._add_cell_scalar(v, k)
            self.extrema[k] = {
                'orig': {'min': str(min(v)), 'max': str(max(v))},
                'user': {'min': str(min(v)), 'max': str(max(v))}
                }
        # supply the combobox with the names to choose from for display
        self.toolbar.cbbx_params.addItems(self.mesh.scalar_names)
        # get the current set parameter
        curr_param = self.toolbar.cbbx_params.currentText()
        # set the first cMin/cMax
        self.toolbar.le_cmin.setText(self.extrema[curr_param]['orig']['min'])
        self.toolbar.le_cmax.setText(self.extrema[curr_param]['orig']['max'])

    def updateParameterView(self, param=None):
        """
        Change the view to given Parameter values.

        Parameter
        ---------
        param: Current text of the just triggered QComboBox

        Note
        ----
        May be overloaded.
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

        if self.toolbar.grp_slice.isChecked():
            x_val = self.toolbar.slice_x.value()
            y_val = self.toolbar.slice_y.value()
            z_val = self.toolbar.slice_z.value()

            # set slicer values into their boxes
            self.toolbar.la_xval.setText(str(round(x_val, 8)))
            self.toolbar.la_yval.setText(str(round(y_val, 8)))
            self.toolbar.la_zval.setText(str(round(z_val, 8)))

            # get the actual slices
            mesh = self.mesh.slice_orthogonal(x=x_val, y=y_val, z=z_val)

        else:
            mesh = self.mesh

        # save the camera position
        # NOTE: this returns [camera position, focal point, and view up]
        self.camera_pos = self.pyvista_widget.camera_position[0]

        # remove the currently displayed mesh
        self.pyvista_widget.remove_actor(self._actor)
        # add the modified one
        self._actor = self.pyvista_widget.add_mesh(mesh, cmap=cMap)

        # update stuff in the toolbar
        self.updateScalarBar()

        # reset the camera position
        self.pyvista_widget.set_position(self.camera_pos)

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
                self.pyvista_widget.update_scalar_bar_range([cmin, cmax], name=param)
        self.pyvista_widget.update()

    def _checkDecimalPoint(self):
        self.toolbar.le_cmin.setText(self.toolbar.le_cmin.text().replace(',', '.'))
        self.toolbar.le_cmax.setText(self.toolbar.le_cmax.text().replace(',', '.'))

    def toggleBbox(self):
        """
        Toggle the visibility of the axis grid surrounding the model.
        """
        checked = not self.toolbar.btn_bbox.isChecked()
        self.pyvista_widget.show_grid(
            show_xaxis=checked,
            show_yaxis=checked,
            show_zaxis=checked,
            show_xlabels=checked,
            show_ylabels=checked,
            show_zlabels=checked
        )
        self.pyvista_widget.update()

    def takeScreenShot(self):
        """
        Save the scene as image.

        Todo
        ----
        + might come in handy to open a dialog where one can choose between 
        black/white background and white/black axis grid and so on
        """
        fname = QFileDialog.getSaveFileName(
            self, 'Open File', None, "Image files (*.jpg *.png)"
            )[0]
        if fname:
            if not len(fname.split('.')) == 2:
                fname += '.png'
            self.pyvista_widget.screenshot(fname)

    def exportMesh(self):
        """
        Save the displayed data as VTK.
        """
        f = QFileDialog.getSaveFileName(
            self, 'Export VTK', None, "VTK file (*.vtk)"
            )[0]
        if f:
            f = f + '.vtk' if not f.lower().endswith('.vtk') else f
            copyfile(self.tmpMesh, f)

    def resetExtrema(self):
        """
        Reset user chosen values to the original ones.
        """
        # get the active scalar/parameter that is displayed currently
        param = self.mesh.active_scalar_name
        self.extrema[param]['user']['min'] = self.extrema[param]['orig']['min']
        self.extrema[param]['user']['max'] = self.extrema[param]['orig']['max']

        # display correctly
        self.updateParameterView(param)

    def _enableSlicers(self):
        if self.toolbar.grp_slice.isChecked():
            self.toolbar.slice_x.setEnabled(True)
            self.toolbar.slice_y.setEnabled(True)
            self.toolbar.slice_z.setEnabled(True)
        else:
            self.toolbar.slice_x.setEnabled(False)
            self.toolbar.slice_y.setEnabled(False)
            self.toolbar.slice_z.setEnabled(False)
        self.updateParameterView()

    def _allowSignals(self):
        # connect signals
        self.toolbar.cbbx_params.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.cbbx_cmap.currentTextChanged.connect(self.updateParameterView)
        self.toolbar.btn_reverse.clicked.connect(self.updateParameterView)
        # self.toolbar.btn_plotlog.clicked.connect(self.updateParameterView)
        self.toolbar.btn_bbox.pressed.connect(self.toggleBbox)
        self.toolbar.btn_screenshot.clicked.connect(self.takeScreenShot)
        self.toolbar.btn_exportVTK.clicked.connect(self.exportMesh)
        self.toolbar.btn_apply.clicked.connect(self.updateScalarBar)
        self.toolbar.btn_reset.clicked.connect(self.resetExtrema)
        self.toolbar.grp_slice.clicked.connect(self._enableSlicers)
        self.toolbar.slice_x.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_y.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_z.sliderReleased.connect(self.updateParameterView)
        self.toolbar.le_cmin.editingFinished.connect(self.updateScalarBar)
        self.toolbar.le_cmax.editingFinished.connect(self.updateScalarBar)


class GToolBar(QWidget):
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
        # combobox to choose the given parameter from
        self.cbbx_params = GComboBox("Scroll/Click to select parameter")

        # combobox to choose the colormap
        self.cbbx_cmap = GComboBox("Scroll/Click to select color map")
        for icon, name in self._createPixmap():
            self.cbbx_cmap.addItem(icon, name)
            self.cbbx_cmap.setIconSize(QSize(40, 15))

        # checkable button to reverse the chosen color scheme
        self.btn_reverse = GButton(
            text="Reverse",
            tooltip="Reverse the chosen color map",
            checkable=True,
            # size=[24, 24]
        )

        # # button for logarithmic values
        # self.btn_plotlog = GButton(
        #     text="Logarithmic",
        #     tooltip="Take logarithmic values of the currently displayed parameter",
        #     checkable=True,
        # )

        # button to make bounding box visible
        self.btn_bbox = GButton(
            text="Toggle Grid",
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

        # button to take a screenshot
        self.btn_screenshot = GButton(
            text="Screenshot",
            tooltip="Save screenshot of the current scene"
            )

        self.btn_exportVTK = GButton(
            text="Export VTK",
            tooltip="Save displayed mesh as VTK"
            )

        # parameter choosing
        lyt_v1 = QVBoxLayout()
        lyt_v1.addWidget(self.cbbx_params)
        lyt_v1.setContentsMargins(2, 2, 2, 2)
        grp_param = QGroupBox("Parameter")
        grp_param.setLayout(lyt_v1)

        # colormap
        lyt_h1 = QVBoxLayout()
        lyt_h1.addWidget(self.cbbx_cmap)
        lyt_h1.addWidget(self.btn_reverse)
        # lyt_h1.addWidget(self.btn_plotlog)
        lyt_h1.setContentsMargins(2, 2, 2, 2)
        grp_cmap = QGroupBox("Color Map")
        grp_cmap.setLayout(lyt_h1)

        # limits
        lyt_h3 = QHBoxLayout()
        lyt_h3.addWidget(QLabel('Min'))
        lyt_h3.addWidget(self.le_cmin)
        lyt_h4 = QHBoxLayout()
        lyt_h4.addWidget(QLabel('Max'))
        lyt_h4.addWidget(self.le_cmax)
        lyt_h5 = QHBoxLayout()
        lyt_h5.addWidget(self.btn_apply)
        lyt_h5.addWidget(self.btn_reset)
        lyt_v3 = QVBoxLayout()
        lyt_v3.addLayout(lyt_h3)
        lyt_v3.addLayout(lyt_h4)
        lyt_v3.addLayout(lyt_h5)
        lyt_v3.setContentsMargins(2, 2, 2, 2)
        grp_limits = QGroupBox("Limits")
        grp_limits.setLayout(lyt_v3)

        # slicer
        lyt_v2 = QVBoxLayout()
        hx = QHBoxLayout()
        hx.addWidget(QLabel("x:"))
        hx.addWidget(self.slice_x)
        self.la_xval = QLineEdit()
        self.la_xval.setReadOnly(True)
        self.la_xval.setFixedWidth(100)
        hx.addWidget(self.la_xval)
        hy = QHBoxLayout()
        hy.addWidget(QLabel("y:"))
        hy.addWidget(self.slice_y)
        self.la_yval = QLineEdit()
        self.la_yval.setReadOnly(True)
        self.la_yval.setFixedWidth(100)
        hy.addWidget(self.la_yval)
        hz = QHBoxLayout()
        hz.addWidget(QLabel("z:"))
        hz.addWidget(self.slice_z)
        self.la_zval = QLineEdit()
        self.la_zval.setReadOnly(True)
        self.la_zval.setFixedWidth(100)
        hz.addWidget(self.la_zval)
        lyt_v2.addLayout(hx)
        lyt_v2.addLayout(hy)
        lyt_v2.addLayout(hz)
        lyt_v2.setContentsMargins(2, 2, 2, 2)
        self.grp_slice = QGroupBox("Slicing")
        self.grp_slice.setCheckable(True)
        self.grp_slice.setChecked(False)
        self.grp_slice.setLayout(lyt_v2)

        # export area
        lyt_h2 = QHBoxLayout()
        lyt_h2.addWidget(self.btn_screenshot)
        lyt_h2.addWidget(self.btn_exportVTK)
        lyt_h2.setContentsMargins(2, 2, 2, 2)
        grp_export = QGroupBox("Export")
        grp_export.setLayout(lyt_h2)

        # widget for the toolbar to better organize the widgets
        lt = QVBoxLayout()
        lt.addWidget(grp_param)
        lt.addWidget(grp_cmap)
        lt.addWidget(grp_limits)
        lt.addWidget(self.grp_slice)
        lt.addStretch(1)
        lt.addWidget(self.btn_bbox)
        lt.addWidget(grp_export)
        lt.setContentsMargins(0, 0, 0, 0)

        self.setLayout(lt)

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
        self.setFixedWidth(200)
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


if __name__ == '__main__':
    pass
