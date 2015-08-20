'''
Created on Jul 27, 2015

@author: hsorby
'''
import re
import subprocess

from glob import glob

from opencmiss.zinc.context import Context
from opencmiss.zinc.streamregion import StreaminformationRegion

from mapclientplugins.monodomain2dstep.utils import sort_numerical_order
import os.path
from opencmiss.zinc.field import Field
from opencmiss.zinc.glyph import Glyph
from opencmiss.zinc.status import OK

rx = re.compile('Time_2_([0-9]+).part0.exnode')

class Mono2DModel(object):
    '''
    classdocs
    '''


    def __init__(self, context_name):
        '''
        Constructor
        '''
        self._context = Context(context_name)
        self._timeKeeper = self._context.getTimekeepermodule().getDefaultTimekeeper()
        glyphmodule = self._context.getGlyphmodule()
        glyphmodule.defineStandardGlyphs()
        materialmodule = self._context.getMaterialmodule()
        materialmodule.defineStandardMaterials()
        tessellationmodule = self._context.getTessellationmodule()
        default_tessellation = tessellationmodule.getDefaultTessellation()
        default_tessellation.setMinimumDivisions(12)
        self._setupSpectrum()
        self.clear()
        
    def clear(self):
        self._min_time = 0.0
        self._max_time = 300.0
        self._step_size = 0.1
        self._time_step = 1.0
        self._x_dis = 4
        self._y_dis = 4
        
        self.clearVisualisation()
        self._region = None
        self._location = None
        self._iron_path = None
        
        
    def initialise(self):
        self._region = self._context.createRegion()
    
    def getDis(self):
        return [self._x_dis, self._y_dis]
        
    def getTimeStep(self):
        return self._time_step
    
    def getStepSize(self):
        return self._step_size
    
    def getMinTime(self):
        return self._min_time
    
    def getMaxTime(self):
        return self._max_time
        
    def getContext(self):
        return self._context
    
    def getRegion(self):
        return self._region
    
    def setLocation(self, location):
        self._location = location
        
    def setIronPath(self, location):
        self._iron_path = location
    
    def simulate(self, step_size, dis):
        self._step_size = step_size
        self._x_dis = dis[0]
        self._y_dis = dis[1]
        
        proc = subprocess.Popen([os.path.join(self._iron_path, 'iron'), str(step_size), str(dis[0]), str(dis[1])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         proc = subprocess.Popen(['/home/hsorby/work/musculoskeletal-software/test-application/iron', str(step_size), str(dis[0]), str(dis[1])],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         stdout, stderr = proc.communicate()
        proc.communicate()
        self.loadSimulation()
        
    def loadSimulation(self):
        self._readMesh()
        
    def setTime(self, time):
        self._timeKeeper.setTime(time)
    
    def _readMesh(self):
        sir = self._region.createStreaminformationRegion()
        files = glob(os.path.join(self._location, 'Time_2_*.part0.exnode'))
        sort_numerical_order(files)
        sir.createStreamresourceFile(os.path.join(self._location.encode('utf-8'), 'MonodomainExample.part0.exelem'))
        for f in files:
            fr = sir.createStreamresourceFile(f.encode('utf-8'))
            m = rx.search(f)
            if m:
                sir.setResourceAttributeReal(fr, StreaminformationRegion.ATTRIBUTE_TIME, int(m.group(1)))
            else:
                print('Big problem not matching a time!!!!')
                
        self._region.read(sir)
        
    def clearVisualisation(self):
        self._surface = None
        self._lines = None
        
    def createVisualisation(self):
#         print('create visulisation')
        self._surface = self._createSurfaceGraphics()
        self._lines = self._createLineGraphics()
#         fieldmodule = self._region.getFieldmodule()
#         fi = fieldmodule.createFielditerator()
#         f = fi.next()
#         while f.isValid():
#             print(f.getName())
#             f = fi.next()
#         print('created surface graphics', self._surface.isValid())

    def _setupSpectrum(self):
#         scenefiltermodule = scene.getScenefiltermodule()
        spectrummodule = self._context.getSpectrummodule()
        default_spectrum = spectrummodule.getDefaultSpectrum()
#         res, minimum, maximum = scene.getSpectrumDataRange(scenefiltermodule.getDefaultScenefilter(), default_spectrum, 1)
        spectrum_component = default_spectrum.getFirstSpectrumcomponent()
        spectrum_component.setRangeMinimum(-95.0)
        spectrum_component.setRangeMaximum(50.0)
#         spectrum_component.setNumberOfBands(10)
#         spectrum_component.setColourMappingType(spectrum_component.COLOUR_MAPPING_TYPE_BANDED)
        
    def _createSurfaceGraphics(self):
        
        scene = self._region.getScene()
        fieldmodule = self._region.getFieldmodule()
        spectrummodule = self._context.getSpectrummodule()
        default_spectrum = spectrummodule.getDefaultSpectrum()
        # We use the beginChange and endChange to wrap any immediate changes and will
        # streamline the rendering of the scene.
        scene.beginChange()
        
        surface = scene.createGraphicsSurfaces()
        coordinate_field = fieldmodule.findFieldByName('Coordinate')
        vm_field = fieldmodule.findFieldByName('Vm')
        surface.setCoordinateField(coordinate_field)
        surface.setDataField(vm_field)
        surface.setSpectrum(default_spectrum)
        #surface.setExterior(True) # show only exterior surfaces
        # Let the scene render the scene.
        scene.endChange()
        # createSurfaceGraphics end
        return surface
        
    def _createLineGraphics(self):
        scene = self._region.getScene()
        fieldmodule = self._region.getFieldmodule()
        fieldmodule.defineAllFaces()
        materialmodule = scene.getMaterialmodule()
        black = materialmodule.findMaterialByName('black')
        # We use the beginChange and endChange to wrap any immediate changes and will
        # streamline the rendering of the scene.
        scene.beginChange()
        coordinate_field = fieldmodule.findFieldByName('Coordinate')
        lines = scene.createGraphicsLines()
        lines.setCoordinateField(coordinate_field)
        lines.setMaterial(black)
        scene.endChange()
        # createSurfaceGraphics end
        return lines
        
