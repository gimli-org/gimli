"""
BUGS:
+ varying the limits (eg exclude a marker) only has effect on the scalar_bar (color_bar), but not on the display of the vtk
+ there is no good coming back from volumetric slicing atm
"""

from shutil import copyfile
import signal
import sys

import pygimli as pg

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QMainWindow, QFrame, QVBoxLayout, QComboBox, QPushButton,
    QFileDialog, QSplitter, QLabel, QAction, QDialog, QStatusBar
)

from .gwidgets import (
    GToolBar, GButton, GLineEdit, GComboBox, GSlider, GDoubleSpinBox, CMAPS
)
pyvista = pg.optImport('pyvista', requiredFor="proper visualization in 3D")


class Show3D(QMainWindow):

    def __init__(self, tmpMesh, application, parent=None):
        super(Show3D, self).__init__(parent)
        self.tmpMesh = tmpMesh
        # storage for the minima and maxima
        self.extrema = {}
        # setup the menubar
        self.setupMenu()
        self.setupWidget()

        self.application = application

        # signals
        signal.signal(signal.SIGINT, self._signal_handler)
        self.acn_close.triggered.connect(self._signal_handler)
        self.acn_hkeys.triggered.connect(self.showHotKeys)

    def _signal_handler(self, sig, frame=None):
        """
        Stop the GUI on CTRL-C, but not the script it was called from.
        from: https://stackoverflow.com/questions/1112343/how-do-i-capture-sigint-in-python
        """
        sys.stderr.write('\r')
        self.application.exit()

    def setupMenu(self):
        bar = self.menuBar()
        # quit the thing
        self.acn_close = QAction("&Quit", self)
        self.acn_close.setShortcut("Q")
        bar.addAction(self.acn_close)
        # about the viewer and help
        ghelp = bar.addMenu("Help")
        self.acn_hkeys = QAction("Hot Keys")
        ghelp.addAction(self.acn_hkeys)

    def showHotKeys(self):
        d = QDialog()
        textfield = QLabel(
        "q - Close pyGIMLi 3D Viewer\n"
        "v - Isometric camera view\n"
        "w - Switch all datasets to a wireframe representation\n"
        "s - Switch all datasets to a surface representation\n"
        "r - Reset the camera to view all datasets\n"
        "shift+click or middle-click - Pan the rendering scene\n"
        "left click - Rotate the rendering scene in 3D\n"
        "ctrl+click - Rotate the rendering scene in 2D (view-plane)\n"
        "mouse-wheel or right-click - Continuously zoom the rendering scene\n"
        )
        textfield.setFont(QFont('Courier'))
        lyt = QVBoxLayout()
        btn_quit = QPushButton('quit')
        btn_quit.clicked.connect(d.done)
        lyt.addWidget(textfield)
        lyt.addWidget(btn_quit)
        lyt.setContentsMargins(2, 2, 2, 2)
        d.setLayout(lyt)
        d.setWindowTitle("Hot Keys")
        d.exec_()

    def setupWidget(self):
        # create the frame
        self.frame = QFrame()

        # add the pyvista interactor object
        self.pyvista_widget = pyvista.QtInteractor()

        vlayout = QVBoxLayout()
        vlayout.setContentsMargins(0, 0, 0, 0)

        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)

        self.toolbar = GToolBar()

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self.toolbar)
        splitter.addWidget(self.pyvista_widget)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 5)

        vlayout.addWidget(splitter)

        self.frame.setLayout(vlayout)

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
        self.statusbar.showMessage("{}".format(self.tmpMesh))
        self.mesh = mesh
        self._actor = self.pyvista_widget.add_mesh(
            mesh, cmap=cMap, show_edges=True)
        self.pyvista_widget.show_bounds(all_edges=True, minor_ticks=True)
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
        pg._d(self.mesh.array_names)
        pg._d([e for e in dir(self.mesh) if 'bar' in e])
        print("-"*60)
        pg._d([e for e in dir(self.pyvista_widget) if 'bar' in e])
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
        This apparently needs to happen when using the gui since on call the
        cell_arrays will be emptied...
        """
        # enable only when there is something to show
        self.toolbar.btn_apply.setEnabled(True)
        self.toolbar.btn_reset.setEnabled(True)
        self.toolbar.spbx_cmin.setEnabled(True)
        self.toolbar.spbx_cmax.setEnabled(True)

        # FIXME: what about the point arrays?!
        for k, v in self.mesh.cell_arrays.items():
            self.mesh._add_cell_array(v, k)
            self.extrema[k] = {
                'orig': {'min': min(v), 'max': max(v)},
                'user': {'min': min(v), 'max': max(v)}
            }
        # supply the combobox with the names to choose from for display
        self.toolbar.cbbx_params.addItems(self.mesh.array_names)
        # get the current set parameter
        curr_param = self.toolbar.cbbx_params.currentText()
        # set the first cMin/cMax
        self.toolbar.spbx_cmin.setValue(
            self.extrema[curr_param]['orig']['min'])
        self.toolbar.spbx_cmax.setValue(
            self.extrema[curr_param]['orig']['max'])

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
            # update the minima and maxima in the limit range
            self.toolbar.spbx_cmin.setRange(
                self.extrema[param]['user']['min'], self.extrema[param]['user']['max'])
            self.toolbar.spbx_cmax.setRange(
                self.extrema[param]['user']['min'], self.extrema[param]['user']['max'])
            self.toolbar.spbx_cmin.setValue(self.extrema[param]['user']['min'])
            self.toolbar.spbx_cmax.setValue(self.extrema[param]['user']['max'])

        cMap = self.toolbar.cbbx_cmap.currentText()
        if self.toolbar.btn_reverse.isChecked():
            cMap += '_r'

        mesh = self.mesh
        if self.toolbar.btn_slice_plane.isChecked() and not self.toolbar.btn_slice_volume.isChecked():
            x_val = self.toolbar.slice_x.value()
            y_val = self.toolbar.slice_y.value()
            z_val = self.toolbar.slice_z.value()

            # set slicer values into their boxes
            self.toolbar.la_xval.setText(str(round(x_val, 8)))
            self.toolbar.la_yval.setText(str(round(y_val, 8)))
            self.toolbar.la_zval.setText(str(round(z_val, 8)))

            # get the actual slices
            mesh = self.mesh.slice_orthogonal(x=x_val, y=y_val, z=z_val)

        if self.toolbar.chbx_threshold.isChecked():
            # get the user defined limits
            cmin = self.toolbar.spbx_cmin.value()
            cmax = self.toolbar.spbx_cmax.value()
            mesh = self.mesh.threshold(value=[cmin, cmax])

        # save the camera position
        # NOTE: this returns [camera position, focal point, and view up]
        self.camera_pos = self.pyvista_widget.camera_position[0]

        # remove the currently displayed mesh
        self.pyvista_widget.remove_actor(self._actor)
        # add the modified one
        if self.toolbar.btn_slice_volume.isChecked() and not self.toolbar.btn_slice_plane.isChecked():
            self._actor = self.pyvista_widget.add_mesh_clip_plane(
                mesh, cmap=cMap, show_edges=True)
        else:
            # in case the plane widget was on.. turn it off
            self.pyvista_widget.disable_plane_widget()
            self._actor = self.pyvista_widget.add_mesh(
            mesh, cmap=cMap, show_edges=True)

        # update stuff in the toolbar
        self.updateScalarBar()

        # reset the camera position
        self.pyvista_widget.set_position(self.camera_pos)

    def updateScalarBar(self):
        """
        When user set limits are made and finished/accepted the color bar
        needs to change.
        """
        cmin = float(self.toolbar.spbx_cmin.value())
        cmax = float(self.toolbar.spbx_cmax.value())
        if cmax >= cmin:
            # get the active scalar/parameter that is displayed currently
            param = self.mesh.active_scalar_name

            # update the user extrema
            self.extrema[param]['user']['min'] = cmin
            self.extrema[param]['user']['max'] = cmax
            # NOTE: has no effect on the displayed vtk
            pg._d("RESET SCALAR BAR LIMITS")
            self.pyvista_widget.update_scalar_bar_range([cmin, cmax])
        self.pyvista_widget.update()

    def toggleBbox(self):
        """
        Toggle the visibility of the axis grid surrounding the model.
        """
        checked = not self.toolbar.btn_bbox.isChecked()
        self.pyvista_widget.show_bounds(
            all_edges=checked,
            minor_ticks=checked,
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
        if self.toolbar.btn_slice_plane.isChecked():
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
        self.toolbar.cbbx_params.currentTextChanged.connect(
            self.updateParameterView)
        self.toolbar.cbbx_cmap.currentTextChanged.connect(
            self.updateParameterView)
        self.toolbar.btn_reverse.clicked.connect(self.updateParameterView)
        # self.toolbar.btn_plotlog.clicked.connect(self.updateParameterView)
        self.toolbar.btn_bbox.pressed.connect(self.toggleBbox)
        self.toolbar.btn_screenshot.clicked.connect(self.takeScreenShot)
        self.toolbar.btn_exportVTK.clicked.connect(self.exportMesh)
        self.toolbar.chbx_threshold.clicked.connect(self.updateParameterView)
        self.toolbar.chbx_threshold.clicked.connect(self._checkStatusThreshold)
        self.toolbar.btn_apply.clicked.connect(self.updateParameterView)
        self.toolbar.btn_reset.clicked.connect(self.resetExtrema)
        self.toolbar.btn_slice_plane.clicked.connect(self._enableSlicers)
        self.toolbar.btn_slice_plane.clicked.connect(self._checkStatusSlice)
        self.toolbar.btn_slice_volume.clicked.connect(self._checkStatusSlice)
        self.toolbar.btn_slice_volume.clicked.connect(self.updateParameterView)
        self.toolbar.slice_x.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_y.sliderReleased.connect(self.updateParameterView)
        self.toolbar.slice_z.sliderReleased.connect(self.updateParameterView)

    def _checkStatusSlice(self):
        """
        Since its either threshold or slice, just disable the other.
        """
        if self.toolbar.btn_slice_plane.isChecked() or self.toolbar.btn_slice_volume.isChecked():
            self.toolbar.chbx_threshold.setChecked(False)
        if self.toolbar.btn_slice_plane.isChecked():
            self.toolbar.btn_slice_volume.setChecked(False)
        elif self.toolbar.btn_slice_volume.isChecked():
            self.toolbar.btn_slice_plane.setChecked(False)

    def _checkStatusThreshold(self):
        """
        Since its either threshold or slice, just disable the other.
        """
        if self.toolbar.chbx_threshold.isChecked():
            self.toolbar.btn_slice_plane.setChecked(False)
            self.toolbar.btn_slice_volume.setChecked(False)


if __name__ == '__main__':
    pass
