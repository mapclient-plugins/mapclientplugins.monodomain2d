'''
Created on Jul 27, 2015

@author: hsorby
'''
import re
import os
import os.path
import sys
import math

from glob import glob


# Intialise OpenCMISS
from opencmiss.iron import iron
from opencmiss.zinc.context import Context
from opencmiss.zinc.streamregion import StreaminformationRegion
from opencmiss.zinc.field import Field
from opencmiss.zinc.glyph import Glyph
from opencmiss.zinc.status import OK

from mapclientplugins.monodomain2dstep.utils import sort_numerical_order


rx = re.compile('Time_2_([0-9]+).part0.exnode')


class Mono2DModel(object):


    def __init__(self, context_name):
        '''
        Constructor
        '''
        self._simulation_root = None
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
        self._time_step = 0.1
        self._x_dis = 4
        self._y_dis = 4
        
        self.clearVisualisation()
        self._region = None
        self._location = None
        self._iron_path = None
        
        
    def initialise(self):
        del self._region
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

    def setSimulationRoot(self, location):
        self._simulation_root = location
    
    def simulate(self, step_size, dis):
        self._step_size = step_size
        self._x_dis = dis[0]
        self._y_dis = dis[1]
        cellml_file = os.path.join(self._simulation_root, 'cellml-models', 'n98.xml')

        # pde_time_step, end_time, output_freq, cellml_file, #x el, #y el
        cwd = os.getcwd()
        os.chdir(self._location)
        solve(0.01, 3.0, 10, cellml_file, self._x_dis, self._y_dis)
        os.chdir(cwd)
        
    def loadSimulation(self):
        self._readMesh()
        
    def setTime(self, time):
        self._timeKeeper.setTime(time)
    
    def _readMesh(self):
        sir = self._region.createStreaminformationRegion()
        files = glob(os.path.join(self._location, 'Time_2_*.part0.exnode'))
        sort_numerical_order(files)
        sir.createStreamresourceFile(os.path.join(self._location, 'MonodomainExample.part0.exelem'))
        for f in files:
            fr = sir.createStreamresourceFile(f)
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
        # fieldmodule = self._region.getFieldmodule()
        # fi = fieldmodule.createFielditerator()
        # f = fi.next()
        # while f.isValid():
        #     print(f.getName())
        #     f = fi.next()
        # print('created surface graphics', self._surface.isValid())

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

def solve(pde_time_step, end_time, output_freq, cellml_file, num_x_els, num_y_els):
    # Set problem parameters
    #DOC-START parameters
    # 2D domain size
    height = 1.0
    width = 1.0
    numberOfXElements = num_x_els
    numberOfYElements = num_y_els

    # Materials parameters
    Am = 193.6
    Cm = 0.014651
    conductivity = 0.1

    # Simulation parameters
    stimValue = 100.0
    stimStop = 0.1
    timeStop = end_time
    odeTimeStep = 0.00001
    pdeTimeStep = pde_time_step
    outputFrequency = output_freq
    #DOC-END parameters

    #Setup field number handles
    coordinateSystemUserNumber = 1
    regionUserNumber = 2
    basisUserNumber = 1
    pressureBasisUserNumber = 2
    generatedMeshUserNumber = 1
    meshUserNumber = 1
    cellMLUserNumber = 1
    decompositionUserNumber = 1
    equationsSetUserNumber = 1
    problemUserNumber = 1
    #Mesh component numbers
    linearMeshComponentNumber = 1
    #Fields
    geometricFieldUserNumber = 1
    fibreFieldUserNumber = 2
    dependentFieldUserNumber = 3
    materialsFieldUserNumber = 4
    equationsSetFieldUserNumber = 5
    cellMLModelsFieldUserNumber = 6
    cellMLStateFieldUserNumber = 7
    cellMLParametersFieldUserNumber = 8
    cellMLIntermediateFieldUserNumber = 9

    #DOC-START parallel information
    # Get the number of computational nodes and this computational node number
    numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
    computationalNodeNumber = iron.ComputationalNodeNumberGet()
    #DOC-END parallel information

    #DOC-START initialisation
    # Create a 2D rectangular cartesian coordinate system
    coordinateSystem = iron.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.DimensionSet(2)
    coordinateSystem.CreateFinish()

    # Create a region and assign the coordinate system to the region
    region = iron.Region()
    region.CreateStart(regionUserNumber,iron.WorldRegion)
    region.LabelSet("Region")
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()
    #DOC-END initialisation

    #DOC-START basis
    # Define a bilinear Lagrange basis
    basis = iron.Basis()
    basis.CreateStart(basisUserNumber)
    basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    basis.numberOfXi = 2
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
    basis.quadratureNumberOfGaussXi = [3]*2
    basis.CreateFinish()
    #DOC-END basis

    #DOC-START generated mesh
    # Create a generated mesh
    generatedMesh = iron.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber,region)
    generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
    generatedMesh.basis = [basis]
    generatedMesh.extent = [width,height]
    generatedMesh.numberOfElements = [numberOfXElements,numberOfYElements]

    mesh = iron.Mesh()
    generatedMesh.CreateFinish(meshUserNumber,mesh)
    #DOC-END generated mesh

    #DOC-START decomposition
    # Create a decomposition for the mesh
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decompositionUserNumber,mesh)
    decomposition.type = iron.DecompositionTypes.CALCULATED
    decomposition.numberOfDomains = numberOfComputationalNodes
    decomposition.CreateFinish()
    #DOC-END decomposition

    #DOC-START geometry
    # Create a field for the geometry
    geometricField = iron.Field()
    geometricField.CreateStart(geometricFieldUserNumber, region)
    geometricField.meshDecomposition = decomposition
    geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Coordinate")
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, linearMeshComponentNumber)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 2, linearMeshComponentNumber)
    geometricField.CreateFinish()

    # Set geometry from the generated mesh
    generatedMesh.GeometricParametersCalculate(geometricField)
    #DOC-END geometry

    #DOC-START equations set
    # Create the equations_set
    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    equationsSetSpecification = [iron.EquationsSetClasses.BIOELECTRICS,
            iron.EquationsSetTypes.MONODOMAIN_EQUATION,
            iron.EquationsSetSubtypes.NONE]
    equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
            equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()
    #DOC-END equations set

    #DOC-START equations set fields
    # Create the dependent Field
    dependentField = iron.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
    equationsSet.DependentCreateFinish()

    # Create the materials Field
    materialsField = iron.Field()
    equationsSet.MaterialsCreateStart(materialsFieldUserNumber, materialsField)
    equationsSet.MaterialsCreateFinish()

    # Set the materials values
    # Set Am
    materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,Am)
    # Set Cm
    materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,Cm)
    # Set conductivity
    materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,conductivity)
    materialsField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,conductivity)
    #DOC-END equations set fields

    #DOC-START create cellml environment
    # Create the CellML environment
    cellML = iron.CellML()
    cellML.CreateStart(cellMLUserNumber, region)
    # Import a Nobel 98 cell model from a file
    noble98Model = cellML.ModelImport(cellml_file)
    #DOC-END create cellml environment

    #DOC-START flag variables
    # Now we have imported the model we are able to specify which variables from the model we want to set from openCMISS
    cellML.VariableSetAsKnown(noble98Model, "fast_sodium_current/g_Na")
    cellML.VariableSetAsKnown(noble98Model, "membrane/IStim")
    # and variables to get from the CellML 
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_K1")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_to")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_K")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_K_ATP")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_Ca_L_K")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_b_K")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_NaK")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_Na")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_b_Na")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_Ca_L_Na")
    cellML.VariableSetAsWanted(noble98Model, "membrane/i_NaCa")
    #DOC-END flag variables

    #DOC-START create cellml finish
    cellML.CreateFinish()
    #DOC-END create cellml finish

    #DOC-START map Vm components
    # Start the creation of CellML <--> OpenCMISS field maps
    cellML.FieldMapsCreateStart()
    #Now we can set up the field variable component <--> CellML model variable mappings.
    #Map Vm
    cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U,1, iron.FieldParameterSetTypes.VALUES,noble98Model,"membrane/V", iron.FieldParameterSetTypes.VALUES)
    cellML.CreateCellMLToFieldMap(noble98Model,"membrane/V", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)

    #Finish the creation of CellML <--> OpenCMISS field maps
    cellML.FieldMapsCreateFinish()

    # Set the initial Vm values
    dependentField.ComponentValuesInitialise(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,-92.5)
    #DOC-END map Vm components

    #DOC-START define CellML models field
    #Create the CellML models field
    cellMLModelsField = iron.Field()
    cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber, cellMLModelsField)
    cellML.ModelsFieldCreateFinish()
    #DOC-END define CellML models field

    #DOC-START define CellML state field
    #Create the CellML state field 
    cellMLStateField = iron.Field()
    cellML.StateFieldCreateStart(cellMLStateFieldUserNumber, cellMLStateField)
    cellML.StateFieldCreateFinish()
    #DOC-END define CellML state field

    #DOC-START define CellML parameters and intermediate fields
    #Create the CellML parameters field 
    cellMLParametersField = iron.Field()
    cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber, cellMLParametersField)
    cellML.ParametersFieldCreateFinish()

    #  Create the CellML intermediate field 
    cellMLIntermediateField = iron.Field()
    cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber, cellMLIntermediateField)
    cellML.IntermediateFieldCreateFinish()
    #DOC-END define CellML parameters and intermediate fields

    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Find the domains of the first and last nodes
    firstNodeNumber = 1
    lastNodeNumber = (numberOfXElements+1)*(numberOfYElements+1)
    firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber, 1)
    lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber, 1)

    # Set the stimulus on half the bottom nodes
    stimComponent = cellML.FieldComponentGet(noble98Model, iron.CellMLFieldTypes.PARAMETERS, "membrane/IStim")
    for node in range(1,int(numberOfXElements/2.0 + 0.5) + 1):
        nodeDomain = decomposition.NodeDomainGet(node,1)
        if nodeDomain == computationalNodeNumber:
            cellMLParametersField.ParameterSetUpdateNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, stimComponent, stimValue)

    # Set up the gNa gradient
    gNaComponent = cellML.FieldComponentGet(noble98Model, iron.CellMLFieldTypes.PARAMETERS, "fast_sodium_current/g_Na")
    for node in range(1,lastNodeNumber):
        nodeDomain = decomposition.NodeDomainGet(node,1)
        if nodeDomain == computationalNodeNumber:
            x = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, 1)
            y = geometricField.ParameterSetGetNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, 2)
            distance = math.sqrt(x*x + y*y)/math.sqrt(width*width + height*height)
            gNaValue = 2*(distance + 0.5)*0.3855
            cellMLParametersField.ParameterSetUpdateNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, gNaComponent, gNaValue)

    #DOC-START define monodomain problem
    #Define the problem
    problem = iron.Problem()
    problemSpecification = [iron.ProblemClasses.BIOELECTRICS,
        iron.ProblemTypes.MONODOMAIN_EQUATION,
        iron.ProblemSubtypes.MONODOMAIN_GUDUNOV_SPLIT]
    problem.CreateStart(problemUserNumber, problemSpecification)
    problem.CreateFinish()
    #DOC-END define monodomain problem

    #Create the problem control loop
    problem.ControlLoopCreateStart()
    controlLoop = iron.ControlLoop()
    problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
    controlLoop.TimesSet(0.0,stimStop,pdeTimeStep)
    controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.TIMING)
    controlLoop.TimeOutputSet(outputFrequency)
    problem.ControlLoopCreateFinish()

    #Create the problem solvers
    daeSolver = iron.Solver()
    dynamicSolver = iron.Solver()
    problem.SolversCreateStart()
    # Get the first DAE solver
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,daeSolver)
    daeSolver.DAETimeStepSet(odeTimeStep)
    daeSolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
    # Get the second dynamic solver for the parabolic problem
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,dynamicSolver)
    dynamicSolver.OutputTypeSet(iron.SolverOutputTypes.NONE)
    problem.SolversCreateFinish()

    #DOC-START define CellML solver
    #Create the problem solver CellML equations
    cellMLEquations = iron.CellMLEquations()
    problem.CellMLEquationsCreateStart()
    daeSolver.CellMLEquationsGet(cellMLEquations)
    cellmlIndex = cellMLEquations.CellMLAdd(cellML)
    problem.CellMLEquationsCreateFinish()
    #DOC-END define CellML solver

    #Create the problem solver PDE equations
    solverEquations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    dynamicSolver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Prescribe any boundary conditions 
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
    solverEquations.BoundaryConditionsCreateFinish()

    # Solve the problem until stimStop
    problem.Solve()

    # Now turn the stimulus off
    for node in range(1,int(numberOfXElements/2.0 + 0.5) + 1):
        nodeDomain = decomposition.NodeDomainGet(node,1)
        if nodeDomain == computationalNodeNumber:
            cellMLParametersField.ParameterSetUpdateNode(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1, node, stimComponent, 0.0)

    #Set the time loop from stimStop to timeStop
    controlLoop.TimesSet(stimStop,timeStop,pdeTimeStep)

    # Now solve the problem from stim stop until time stop
    problem.Solve()

    # Export the results, here we export them as standard exnode, exelem files
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("MonodomainExample","FORTRAN")
    fields.ElementsExport("MonodomainExample","FORTRAN")
    fields.Finalise()

    coordinateSystem.Destroy()
    region.Destroy()
    basis.Destroy()
    problem.Destroy()
    
