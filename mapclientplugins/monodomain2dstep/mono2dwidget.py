'''
Created on Jul 27, 2015

@author: hsorby
'''
import os.path

from PySide import QtGui, QtCore

from mapclientplugins.monodomain2dstep.ui_mono2dwidget import Ui_Mono2DWidget

from mapclientplugins.monodomain2dstep.mono2dmodel import Mono2DModel

class Mono2DWidget(QtGui.QWidget):
    '''
    classdocs
    '''

    def __init__(self, parent=None):
        '''
        Constructor
        '''
        super(Mono2DWidget, self).__init__(parent)
        self._ui = Ui_Mono2DWidget()
        self._ui.setupUi(self)
        self._ui.doubleSpinBoxStepSize.setVisible(False)
        self._ui.labelStepSize.setVisible(False)
        # The number of decimals defines the minimum positive
        # number we can set.
        self._ui.doubleSpinBoxStepSize.setDecimals(8)
        self._ui.doubleSpinBoxStepSize.setMinimum(0.00000001)
        
        self._model_converged = Mono2DModel('converged')
        self._model_experimental = Mono2DModel('experimental')
        
        self._ui.widgetSceneviewerConverged.setContext(self._model_converged.getContext())
        self._ui.widgetSceneviewerExperimental.setContext(self._model_experimental.getContext())
        
        self._callback = None
        self._timer = QtCore.QTimer()
        
        self._makeConnections()
        
    def clear(self):
        self._model_converged.clear()
        self._model_experimental.clear()
        
    def initialise(self, data_location):
        self._model_converged.initialise()
        self._model_experimental.initialise()
        self._ui.pushButtonPlayStop.setText('Play')
        self._setSliderValues()
        self._ui.labelTime.setText('{:10.4f}'.format(self._model_experimental.getMinTime()))
        self._ui.doubleSpinBoxStepSize.setValue(self._model_experimental.getStepSize())
        dis = self._model_experimental.getDis()
        self._ui.spinBoxXDiscretisation.setValue(dis[0])
        self._ui.spinBoxYDiscretisation.setValue(dis[1])
        self._timer.setInterval(0) #*self._model_experimental.getTimeStep())

        
        self._model_converged.setLocation(os.path.join(data_location, 'Monodomain2D/converged'))
        self._model_experimental.setLocation(os.path.join(data_location, 'Monodomain2D/experimental'))
        self._model_converged.setIronPath(os.path.join(data_location, 'bin'))
        self._model_experimental.setIronPath(os.path.join(data_location, 'bin'))

        self._model_converged.simulate(0.1, [7, 7])
        self._model_converged.createVisualisation()
#         self._model_experimental.simulate(0.1, [dis[0], dis[1]])
#         self._model_experimental.createVisualisation()

    def _tweakView(self):
        p = self._ui.widgetSceneviewerConverged.getViewParameters()
        print(p)

    def _makeConnections(self):
        self._ui.pushButtonContinue.clicked.connect(self._continueClicked)
        self._ui.pushButtonSimulate.clicked.connect(self._simulateClicked)
        self._ui.pushButtonPlayStop.clicked.connect(self._playStopClicked)
        self._timer.timeout.connect(self._timerTimedOut)
        self._ui.horizontalSliderTime.valueChanged.connect(self._timeChanged)
        self._ui.widgetSceneviewerConverged.graphicsInitialized.connect(self._graphicsInitialized)
        self._ui.widgetSceneviewerExperimental.graphicsInitialized.connect(self._graphicsInitialized)
        
    def _graphicsInitialized(self):
        sender = self.sender()
        sceneviewer = sender.getSceneviewer()
        if sceneviewer is not None:
            if sender is self._ui.widgetSceneviewerConverged:
                scene = self._model_converged.getRegion().getScene()
            elif sender is self._ui.widgetSceneviewerExperimental:
                scene = self._model_experimental.getRegion().getScene()
            sceneviewer.setScene(scene)
            sceneviewer.viewAll()
            sceneviewer.setPerturbLinesFlag(True)
            # We need to tweak the view slightly so that we
            # can see the lines of the elements.
            _, v = sceneviewer.getEyePosition()
            v[1] += 0.01
            sceneviewer.setEyePosition(v)
        
    def _continueClicked(self):
        self._callback()
        
    def registerCallback(self, callback):
        self._callback = callback
        
    def _simulateClicked(self):
        x_dis = self._ui.spinBoxXDiscretisation.value()
        y_dis = self._ui.spinBoxYDiscretisation.value()
        step_size = self._ui.doubleSpinBoxStepSize.value()
        
        QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
#         self._model_experimental.clearVisualisation()
        self._model_experimental.simulate(step_size, [x_dis, y_dis])
        self._model_experimental.createVisualisation()
        self._ui.widgetSceneviewerExperimental.viewAll()
        QtGui.QApplication.restoreOverrideCursor()
        
    def _setSliderValues(self):
        step_size = self._model_experimental.getStepSize()
        slider_range = (self._model_experimental.getMaxTime() - self._model_experimental.getMinTime())/step_size
        self._ui.horizontalSliderTime.setRange(0, slider_range)
        self._ui.horizontalSliderTime.setValue(0)
        
    def _setTime(self, value):
        """
        The value here is the slider value and not the actual desired time
        we need to convert before using it.
        """
        step_size = self._model_experimental.getStepSize()
        time = self._model_experimental.getMinTime() + value * step_size
        self._ui.labelTime.setText('{:10.4f}'.format(time))
        self._model_converged.setTime(10.0*time)
        self._model_experimental.setTime(time)
        
    def _timeChanged(self, value):
        """
        Deals with events from the user manually changing the slider 
        to a new value.
        """
        self._setTime(value)
        
    def _timerTimedOut(self):
        """
        Deals with timeout events triggered by the timer. i.e. the 
        Play button is active.
        """
        value = self._ui.horizontalSliderTime.value()
        max_value = self._ui.horizontalSliderTime.maximum()
        value += 1
        if max_value < value:
            value = 0

        self._ui.horizontalSliderTime.setValue(value)
        self._setTime(value)
        
    def _playStopClicked(self):
        button_text = self.sender().text()
        if button_text == 'Play':
            self._timer.start()
            self._ui.pushButtonPlayStop.setText('Stop')
            self._ui.horizontalSliderTime.setEnabled(False)
        else:
            self._timer.stop()
            self._ui.pushButtonPlayStop.setText('Play')
            self._ui.horizontalSliderTime.setEnabled(True)
        
        
        