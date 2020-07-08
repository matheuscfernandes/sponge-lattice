from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np
import math

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

execfile('Functions.py')

def RunSimulation(LengthOfBridge,THICKNESS):
    Mdb()

    PERIMITERRATIO=3.25
    spacing = 1.
    TOL = 10E-7
    DD = 0.1*spacing
    YoungsModulus = 1.
    PoissonsRatio = 0.3
    AppYDisp = -2*spacing
    x2 = spacing/(sqrt(2)+2)
    Lambda = sqrt(2.)

    Model_DimpleStr=mdb.models['Model-1']

    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    #CREATING THE MAIN UNICEL
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        2*spacing, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 2*spacing), point2=(
        0, 0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
        spacing, 2*spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
        2*spacing, spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 2*spacing), point2=(
        2*spacing, 2*spacing))
    # CREATING UNICEL FOR THE ONE DIAGONAL PART
    if DesignB:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing, 0), point2=(
            0, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
            2*spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 2*spacing))

    if DesignC:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing, 0), point2=(
            0, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
            2*spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 2*spacing))

        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 0), point2=(
            2*spacing, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 0))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 2*spacing), point2=(
            2*spacing, spacing))

    # CREATING UNICEL FOR THE TWO DIAGONAL PART
    if DesignA:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing-x2, 0), point2=(
            2*spacing, spacing+x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing+x2, 0), point2=(
            2*spacing, spacing-x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing-x2, 0), point2=(
            0, 2*spacing-x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing, x2), point2=(
            x2, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing-x2), point2=(
            spacing+x2, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing+x2), point2=(
            spacing-x2, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, x2), point2=(
            x2, 0))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing-x2, 2*spacing), point2=(
            2*spacing, 2*spacing-x2))

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseWire(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    ## CREATING SINGLE FLOATING CELL FOR THE LAST ODD NUMBER ENDING

    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    #CREATING THE MAIN UNICEL
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        spacing, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 2*spacing), point2=(
        0, 0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
        spacing, 2*spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
        spacing, spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 2*spacing), point2=(
        spacing, 2*spacing))
    # CREATING UNICEL FOR THE ONE DIAGONAL PART
    if DesignB:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, spacing), point2=(
            0, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 2*spacing))

    if DesignC:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 2*spacing), point2=(
            spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 2*spacing))

        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 0), point2=(
            spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
            spacing, 0))
        

    # CREATING UNICEL FOR THE TWO DIAGONAL PART
    if DesignA:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing-x2, 0), point2=(
            spacing, x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, spacing-x2), point2=(
            0, 2*spacing-x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, spacing+x2), point2=(
            x2, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing-x2), point2=(
            spacing, 2*spacing-x2))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing+x2), point2=(
            spacing-x2, 2*spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, x2), point2=(
            x2, 0))

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-2', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-2'].BaseWire(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']


    #CREATING SETS FOR DIAGONAL CROSSECTION AND OTHER CROSSECTIONS
    Part_Full = mdb.models['Model-1'].parts['Part-1']
    edgeIndex=[]
    middleIndex=[]
    diagIndex=[]
    for ii in range(len(Part_Full.edges)):
        pointOnX=Part_Full.edges[ii].pointOn[0][0]
        pointOnY=Part_Full.edges[ii].pointOn[0][1]
        if pointOnX==spacing or pointOnY==spacing:
            middleIndex.append(Part_Full.edges[ii:ii+1])
        elif pointOnX==0.0 or pointOnX==2*spacing or pointOnY==0.0 or pointOnY==2*spacing:
            edgeIndex.append(Part_Full.edges[ii:ii+1])
        else:
            diagIndex.append(Part_Full.edges[ii:ii+1])
    edgeIndex=tuple(edgeIndex)
    middleIndex=tuple(middleIndex)
    diagIndex=tuple(diagIndex)

    if DesignA or DesignB or DesignC: Part_Full.Set(name='DIAGONAL', edges=diagIndex)
    Part_Full.Set(name='EDGES', edges=edgeIndex+middleIndex)
    # Part_Full.Set(name='MIDDLE', edges=middleIndex)

    #FOR THE OUTSIDE EDGES
    Part_Full = mdb.models['Model-1'].parts['Part-2']
    edgeIndex=[]
    middleIndex=[]
    diagIndex=[]
    for ii in range(len(Part_Full.edges)):
        pointOnX=Part_Full.edges[ii].pointOn[0][0]
        pointOnY=Part_Full.edges[ii].pointOn[0][1]
        if pointOnX==spacing or pointOnY==spacing:
            middleIndex.append(Part_Full.edges[ii:ii+1])
        elif pointOnX==0.0 or pointOnX==2*spacing or pointOnY==0.0 or pointOnY==2*spacing:
            edgeIndex.append(Part_Full.edges[ii:ii+1])
        else:
            diagIndex.append(Part_Full.edges[ii:ii+1])
    edgeIndex=tuple(edgeIndex)
    middleIndex=tuple(middleIndex)
    diagIndex=tuple(diagIndex)

    if DesignA or DesignB or DesignC: Part_Full.Set(name='DIAGONAL', edges=diagIndex)
    Part_Full.Set(name='EDGES', edges=edgeIndex+middleIndex)
    # Part_Full.Set(name='MIDDLE', edges=middleIndex)
    

    #CREATING INSTANCE AND MESHING
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    Instant_Full=mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    

    #CREATING LINEAR PATTERN FOR FINITE SIZE
    mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
        0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Part-1-1', ), number1=LengthOfBridge, 
        number2=1, spacing1=2.*spacing, spacing2=0)


    
    #ADD OUTER EDGE PART 2 INTO THE ASSEMBLY
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=mdb.models['Model-1'].parts['Part-2'])

    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-1', ), 
    vector=(2*spacing*LengthOfBridge, 0.0, 0.0))

    mdb.models['Model-1'].rootAssembly.rotate(angle=180.0, axisDirection=(0.0, 0.0, 
    1.0), axisPoint=(spacing*LengthOfBridge+spacing/2., spacing, 0.0), instanceList=tuple(mdb.models['Model-1'].rootAssembly.instances.keys()))

    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=BOTH, 
        instances=(mdb.models['Model-1'].rootAssembly.instances.values()), mergeNodes=
        BOUNDARY_ONLY, name='Part-3', nodeMergingTolerance=1e-06, 
        originalInstances=DELETE)
    
    mdb.models['Model-1'].parts.changeKey(fromName='Part-1', toName='Part-1A')
    mdb.models['Model-1'].parts.changeKey(fromName='Part-3', toName='Part-1')
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Part-3-1', 
        toName='Part-1-1')
    
    mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.01, 
        minSizeFactor=0.01, size=0.03)

    #change element type to quadratic interpolation
    mdb.models['Model-1'].parts['Part-1'].Set(name='ALL',edges=mdb.models['Model-1'].parts['Part-1'].edges[:])
    mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
        elemCode=B22, elemLibrary=STANDARD), ), regions=mdb.models['Model-1'].parts['Part-1'].sets['ALL'])

    mdb.models['Model-1'].parts['Part-1'].generateMesh()
    
    #CREATING STEP
    mdb.models['Model-1'].StaticStep(initialInc=0.001, name='Step-1', nlgeom=ON, 
        previous='Initial')

    #CONTACT DEFINITION
    # mdb.models['Model-1'].ContactProperty('IntProp-1')
    # mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    #     allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
    #     pressureOverclosure=HARD)
    # mdb.models['Model-1'].ContactStd(createStepName='Initial', name='Int-1')
    # mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    #     stepName='Initial', useAllstar=ON)
    # mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    #     assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')

    ## CREATING AND WORKING WITH THE VIRTUAL NODES
    #CREATING SETS
    Part_Full = mdb.models['Model-1'].parts['Part-1']

    Part_Full.Set(name='FixedEdges', nodes=
        (Part_Full.nodes.getByBoundingBox(0 - TOL, 0 - TOL, -TOL, + TOL, + TOL, TOL), 
        Part_Full.nodes.getByBoundingBox(LengthOfBridge*2.*spacing +spacing - TOL,  - TOL, -TOL,LengthOfBridge*2.*spacing +spacing+ TOL, + TOL, TOL) )) 


    Part_Full.Set(name='Load', nodes=
        (Part_Full.nodes.getByBoundingBox(LengthOfBridge*spacing +spacing/2.- TOL,  2*spacing- TOL, -TOL,LengthOfBridge*spacing +spacing/2. + TOL, 2*spacing+ TOL, TOL), )) 

   #APPLY BC
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-FIXED', region=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['FixedEdges'], 
        u1=0.0, u2=0.0, ur3=UNSET)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-DISP', region=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['Load'], 
        u1=UNSET, u2=AppYDisp, ur3=UNSET)
   


    # mdb.models['Model-1'].ConcentratedForce(cf1=0.0, cf2=AppYLoad, createStepName=
    #     'Step-1', distributionType=UNIFORM, field='', localCsys=None,  name='BC-LOAD', region=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['Load'])


    # DEFINING MATERIAL PROPERTIES AND SECTION PROPERTIES
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=
        ISOTROPIC, table=((1.0, 0.0), ), testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=POISSON_RATIO)
    # mdb.models['Model-1'].materials['Material-1'].Elastic(table=((1.0, 0.0), ))

    # CREATE PROFILE AND SECTION ASSIGNMENT - THIS IS FOR THE ROUND CROSSECTION
    if DesignB:
        DEDGES = DD
        DDIAGONALS = DEDGES/sqrt(2.)

    elif DesignC:
        DEDGES = DD
        DDIAGONALS = DEDGES/2.

    elif DesignA:
        DEDGES=(Lambda/(2.*(1.+Lambda)))*(1.+1./sqrt(2.))/5.
        DDIAGONALS=(1./(2.*2.*(1.+Lambda)))*(1.+sqrt(2.))/5.

    else:
        DEDGES = (DD*(1.+sqrt(2.)/2.))
    
    #DEFINING SECTION FOR THE GLOBAL EDGES
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES, name='EDGES')

    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='EDGES', poissonRatio=0.0, 
        profile='EDGES', temperatureVar=LINEAR)

    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1'].sets['EDGES'], sectionName=
        'EDGES', thicknessAssignment=FROM_SECTION)

    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1'].sets['EDGES'])

    
    # DEFINING SECTION FOR THE DIAGONAL STRUTS
    if DesignB or DesignA or DesignC:
        mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DDIAGONALS, name='DIAGONAL')
        mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
            DURING_ANALYSIS, material='Material-1', name='DIAGONAL', poissonRatio=0.0, 
            profile='DIAGONAL', temperatureVar=LINEAR)
        mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=
            mdb.models['Model-1'].parts['Part-1'].sets['DIAGONAL'], sectionName=
            'DIAGONAL', thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=
            N1_COSINES, n1=(0.0, 0.0, -1.0), region=
            mdb.models['Model-1'].parts['Part-1'].sets['DIAGONAL'])

    # CREATE JOB AND RUN IT
    if DesignB:
        JobName="DesignB"
    elif DesignA:
        JobName="DesignA"
    elif DesignC:
        JobName="DesignC"
    else:
        JobName="DesignD"
    
    DeleteAbaqusFiles(JobName)
    mdb.models['Model-1'].rootAssembly.regenerate()
    
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=
        100)
    
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=THREADS, name=JobName, nodalOutputPrecision=SINGLE, 
        numCpus=8, numDomains=8, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs[JobName].submit(consistencyChecking=OFF)
    mdb.jobs[JobName].waitForCompletion()

    ExtractSetRF('PART-1-1','Step-1',JobName,'LOAD')
    DeleteAbaqusFilesButODB(JobName)