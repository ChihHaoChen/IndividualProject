__author__ = 'Chih-Hao'

from PyQt4 import QtCore, QtGui


try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

import sys

from PyCurve_2D_qt import Ui_PyCurve_2D
from Meshimport import mesh
from Interpolation_5RBF import interpolation_compute
import Displacement
import Meshexport
import numpy as np

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import time

'''
To import the python files of previous mesh computation
1. Mesh import
2. Displacement computation
3. Do RBF interpolation on the basis of the computed displacement
4. Mesh export for GMSH reading out
'''


class QMainWindow(QtGui.QMainWindow, Ui_PyCurve_2D):
    timeElapsed = QtCore.pyqtSignal(int)
    doneAlert = QtCore.pyqtSignal(str)
    show_string = QtCore.pyqtSignal(float, float)
    status_alert = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super(QMainWindow, self).__init__(parent)

        self.setupUi(self)
        self.ui_restore()

    @QtCore.pyqtSlot()
    def on_button_clicked(self):
        self.matplotlib_widget1.axis.clear()
        self.matplotlib_widget2.axis.clear()
        self.timer_thread.reset_time()
        self.timer_thread.start()

    @QtCore.pyqtSlot(int)
    def on_myThread_timeElapsed(self, seconds):
        self.ShowTime.display(time.strftime('%H:%M:%S', time.gmtime(seconds)))

    @QtCore.pyqtSlot(str)
    def on_myThread_done(self, text):
        self.timer_thread.terminate()
        if text == "Running!":
            self.ShowState.setStyleSheet('color: red')
        elif text == "Done!":
            self.ShowState.setStyleSheet('color: green')
        elif text == "Idle State":
            self.ShowState.setStyleSheet('color: black')

        self.ShowState.setText(text)

    def ui_restore(self):
        self.setupUi(self)
        self.fontC = QtGui.QFont()
        self.fontC.setFamily(_fromUtf8("Calibri"))
        self.fontC.setPointSize(10)
        self.optimise_option = 0
        self.counter_sec = 0
        self.counter_min = 0
        self.counter_hr = 0
        self.button_group = QtGui.QButtonGroup()
        self.two_dim_filename = ""
        self.one_dim_filename = ""
        self.boundary_inner = ""
        self.boundary_outer = ""
        self.mesh_order = 0
        self.local_RBF_radius = 0
        self.CTPS_radius = 0
        self.save_file_name = ""
        self.counter_sec = 0
        self.counter_min = 0
        self.counter_hr = 0
        self.disp_delta_x = 0
        self.disp_delta_y = 0
        self.disp_pivot_x = 0
        self.disp_pivot_y = 0
        self.disp_rot_degree = 0

        self.timer_thread = timerThread(self)

        self.timer_thread.timeElapsed.connect(self.timeElapsed.emit)
        self.timer_thread.doneAlert.connect(self.doneAlert.emit)
        self.timer_thread.timeElapsed.connect(self.on_myThread_timeElapsed)

        self.connect(self.Browse_1D, QtCore.SIGNAL('clicked()'), self.browse_1dim)
        self.connect(self.Browse_2D, QtCore.SIGNAL('clicked()'), self.browse_2dim)
        self.Generate.clicked.connect(self.on_button_clicked)
        self.Stop.clicked.connect(self.stop_thread)
        self.Reset.clicked.connect(self.reset_thread)

        self.checkedbox()
        self.ShowState.setText("Idle State")
        self.connect(self.Generate, QtCore.SIGNAL('clicked()'), self.setup)
        self.ShowTime.setFont(self.fontC)
        self.ShowTime.display(time.strftime('%H:%M:%S', time.gmtime(0)))
        self.matplotlib_widget1 = matplotlib_widget(self)
        self.matplotlib_widget2 = matplotlib_widget(self)
        self.layout1 = QtGui.QVBoxLayout(self.Graph1)
        self.layout2 = QtGui.QVBoxLayout(self.Graph2)
        self.layout1.addWidget(self.matplotlib_widget1)
        self.layout2.addWidget(self.matplotlib_widget2)

        self.myThread = None

    def setup(self):
        self.boundary_inner = self.InnerBound.text()
        self.boundary_outer = self.OuterBound.text()
        self.checkedbox()
        self.ShowState.setText("Running!")
        self.ShowState.setStyleSheet('color: red')
        self.mesh_order = self.Order.text()
        self.local_RBF_radius = self.Radius.text()
        self.CTPS_radius = self.CTPS_Radius.text()
        self.save_file(self.OutputFile.text())
        self.disp_delta_x = self.delta_x.text()
        self.disp_delta_y = self.delta_y.text()
        self.disp_pivot_x = self.pivot_x.text()
        self.disp_pivot_y = self.pivot_y.text()
        self.disp_rot_degree = self.rot_degree.text()

        self.myThread = myThread(self.boundary_inner, self.boundary_outer, self.mesh_order, self.local_RBF_radius, \
                                self.CTPS_radius, self.two_dim_filename, self.one_dim_filename, self.save_file_name, \
                                self.optimise_option, self.disp_delta_x, self.disp_delta_y, self.disp_pivot_x, \
                                self.disp_pivot_y, self.disp_rot_degree)

        self.myThread.start()

        self.myThread.show_string.connect(self.show_string.emit)
        self.myThread.show_string.connect(self.quality_plot)
        self.myThread.status_alert.connect(self.status_alert.emit)
        self.myThread.status_alert.connect(self.on_myThread_done)

    def stop_thread(self):
        if self.myThread.isRunning():
            self.myThread.terminate()
            self.timer_thread.terminate()

        self.ShowState.setText("Suspended!")
        self.ShowState.setStyleSheet('color: red')

    def reset_thread(self):
        if self.myThread.isRunning():
            self.stop_thread()

        self.setupUi(self)
        self.ui_restore()

    def checkedbox(self):
        self.button_group.addButton(self.CheckMeshQ_skewness)
        self.button_group.addButton(self.CheckMeshQ_nonorth)
        self.button_group.addButton(self.CheckTPS)
        self.button_group.addButton(self.CheckCPC2)
        self.button_group.addButton(self.CheckCPC2b)
        self.button_group.addButton(self.CheckCPC4)
        self.button_group.addButton(self.CheckCPC6)
        self.button_group.setExclusive(True)

        if self.CheckMeshQ_skewness.isChecked():
            self.optimise_option = 1
        elif self.CheckMeshQ_nonorth.isChecked():
            self.optimise_option = 2
        elif self.CheckTPS.isChecked():
            self.optimise_option = 3
        elif self.CheckCPC2.isChecked():
            self.optimise_option = 4
        elif self.CheckCPC2b.isChecked():
            self.optimise_option = 5
        elif self.CheckCPC4.isChecked():
            self.optimise_option = 6
        elif self.CheckCPC6.isChecked():
            self.optimise_option = 7
        else:
            self.optimise_option = 0
            self.ShowState.setText("No Chosen")

    def browse_1dim(self):
        self.one_dim_filename = QtGui.QFileDialog.getOpenFileName(self, "Open Mesh File", "*.msh")
        if self.OneDimFile.findText(self.OneDimFile.currentText()) == -1:
            self.OneDimFile.setFont(self.fontC)
            self.OneDimFile.addItem(_fromUtf8(self.one_dim_filename))
        else:
            self.OneDimFile.clear()
            self.OneDimFile.setFont(self.fontC)
            self.OneDimFile.addItem(_fromUtf8(self.one_dim_filename))

    def browse_2dim(self):
        self.two_dim_filename = QtGui.QFileDialog.getOpenFileName(self, "Open Mesh File", "*.msh")

        if self.TwoDimFile.findText(self.TwoDimFile.currentText()) == -1:
            self.TwoDimFile.setFont(self.fontC)
            self.TwoDimFile.addItem(str(self.two_dim_filename))
        else:
            self.TwoDimFile.clear()
            self.TwoDimFile.setFont(self.fontC)
            self.TwoDimFile.addItem(_fromUtf8(self.two_dim_filename))

    def save_file(self, filename):
        self.save_file_name = filename+".msh"

    @QtCore.pyqtSlot(float, float)
    def quality_plot(self, float1, float2):
        skew_text_ave = float1
        nonorht_text_ave = float2
        meanStrSkew = 'The Average Skewness :'
        meanStrSkew += str("{0:.3f}".format(skew_text_ave))
        meanStrSkew += self.myThread.optimise_method
        meanStrNonOr = 'The Average Non-Orthognality :'
        meanStrNonOr += str("{0:.3f}".format(nonorht_text_ave))
        meanStrNonOr += self.myThread.optimise_method

        self.matplotlib_widget1.axis.grid()
        self.matplotlib_widget1.axis.set_title(meanStrSkew)
        self.matplotlib_widget1.axis.hist(self.myThread.skew_max, bins=50, range=(0, 1))
        self.matplotlib_widget1.figure_canvas.draw()

        self.matplotlib_widget2.axis.grid()
        self.matplotlib_widget2.axis.set_title(meanStrNonOr)
        self.matplotlib_widget2.axis.hist(self.myThread.non_orth_angle, bins=100, range=(0, 90))
        self.matplotlib_widget2.figure_canvas.draw()


class matplotlib_widget(QtGui.QWidget):
    def __init__(self, parent=None):
        super(matplotlib_widget, self).__init__(parent)

        self.figure = Figure()
        self.figure_canvas = FigureCanvas(self.figure)

        self.axis = self.figure.add_subplot(111)
        self.axis.grid()
        self.axis.set_title("")
        self.layout = QtGui.QVBoxLayout(self)
        self.layout.addWidget(self.figure_canvas)


class myThread(QtCore.QThread):
    show_string = QtCore.pyqtSignal(float, float)
    status_alert = QtCore.pyqtSignal(str)

    def __init__(self, boundary_inner, boundary_outer, mesh_order, local_RBF_radius, CTPS_radius, two_dim_filename, \
                 one_dim_filename, save_file_name, optimise_option, delta_x, delta_y, pivot_x, pivot_y, rot_degree, \
                 parent=None):
        super(myThread, self).__init__(parent)
        self.boundary_inner = boundary_inner
        self.boundary_outer = boundary_outer
        self.mesh_order = mesh_order
        self.local_RBF_radius = local_RBF_radius
        self.CTPS_radius = CTPS_radius
        self.two_dim_filename = two_dim_filename
        self.one_dim_filename = one_dim_filename
        self.save_file_name = save_file_name
        self.optimise_option = optimise_option
        self.delta_x = delta_x
        self.delta_y = delta_y
        self.pivot_x = pivot_x
        self.pivot_y = pivot_y
        self.rot_degree = rot_degree

        self.skew_max_ave = 0
        self.non_orth_angle_ave = 0
        self.optimise_method = ""
        self.skew_max = np.array([])
        self.non_orth_angle = np.array([])
        self.optimise_method = ""

    def run(self):
        self.setPriority(QtCore.QThread.LowPriority)
        self.mesh_refinement()
        self.signal_emit()

    def mesh_refinement(self):
        phyLine = self.boundary_inner.split(',')
        phyLineOut = self.boundary_outer.split(',')
        order = float(self.mesh_order)
        dim = 2
        radius = float(self.local_RBF_radius)
        radiusCTPS = float(self.CTPS_radius)
        meshM = mesh(dim, order, phyLine, phyLineOut)
        meshT = mesh(dim, order, phyLine, phyLineOut)
        meshIndex = self.optimise_option
        ##The original 2D mesh
        meshM.openfile(self.two_dim_filename)
        ##The 1D mesh by the generation of curvilinear function in Gmsh
        meshT.openfile(self.one_dim_filename)
        ##The displacment function
        disp = Displacement.displacement(meshM, meshT, float(self.delta_x), float(self.delta_y), float(self.pivot_x), \
                                            float(self.pivot_y), float(self.rot_degree))
        ##The main function of interpolation
        final_mesh = interpolation_compute(meshM, disp, dim, meshIndex, radius, radiusCTPS)

        self.skew_max_ave = final_mesh.skew_max_ave
        self.non_orth_angle_ave = final_mesh.non_orth_angle_ave
        self.skew_max = final_mesh.skew_max
        self.non_orth_angle = final_mesh.nonorth_angle
        self.optimise_method = final_mesh.optimise_method

        Meshexport.meshexport(meshM, final_mesh.interpolation(), self.save_file_name)

    def signal_emit(self):
        self.show_string.emit(self.skew_max_ave, self.non_orth_angle_ave)
        self.status_alert.emit("Done!")


class timerThread(QtCore.QThread):
    timeElapsed = QtCore.pyqtSignal(int)
    doneAlert = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super(timerThread, self).__init__(parent)
        self.timeStart = None
        self.time_finish = 0

    def run(self):
        self.setPriority(QtCore.QThread.HighestPriority)
        self.reset_time()
        while True:
            time.sleep(1)
            self.timeElapsed.emit(time.time() - self.timeStart)
            self.doneAlert.emit("Running!")

    def reset_time(self):
        self.timeStart = time.time()


def main():
    app = QtGui.QApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)
    form = QMainWindow()
    form.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()


