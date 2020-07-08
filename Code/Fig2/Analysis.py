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




# CREATE JOB AND RUN IT
if DesignB:
    JobName="DesignB"
    pert=[0.01]
    modes=[1]

elif DesignA:
    if Optimization:
        JobName="DesignA_Optimization"
        pert=[0.01]
        modes=[1]
    else:
        JobName="DesignA"
        pert=[0.01]
        modes=[1]
elif DesignC:
    JobName="DesignC"
    pert=[0.01]
    modes=[1]
else:
    JobName="DesignD"
    pert=[0.01]
    modes=[1]

def RunBuckling(numberofunicells,THICKNESS):
    Mdb()

    PERIMITERRATIO=2
    spacing = 1.
    TOL = 10E-7
    DD = 0.1*spacing
    YDisplacement = -0.01*2*spacing*numberofunicells

    Sep=6.88666541e-01
    Lambda=1000*6.84940074e-04

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
        if Optimization: 
            x2 = spacing*Sep
        else: 
            x2=spacing/(sqrt(2)+2)
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

    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(numberofunicells*2.*spacing, 0.0), point2=(
        numberofunicells*2.*spacing, numberofunicells*2.*spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(numberofunicells*2.*spacing, numberofunicells*2.*spacing), point2=
        (0.0, numberofunicells*2.*spacing))
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
    Part_Full.Set(name='EDGES', edges=edgeIndex)
    Part_Full.Set(name='MIDDLE', edges=middleIndex)

    #FOR THE OUTSIDE EDGES
    Part_Full = mdb.models['Model-1'].parts['Part-2']
    Part_Full.Set(name='EDGES-PART2', edges=Part_Full.edges[:])

    #CREATING INSTANCE AND MESHING
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    Instant_Full=mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])

    #CREATING LINEAR PATTERN FOR FINITE SIZE
    mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
        0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Part-1-1', ), number1=numberofunicells, 
        number2=numberofunicells, spacing1=2.*spacing, spacing2=2.*spacing)

    #ADD OUTER EDGE PART 2 INTO THE ASSEMBLY
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=mdb.models['Model-1'].parts['Part-2'])

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
    mdb.models['Model-1'].BuckleStep(name='Step-1',
                                   numEigen=6, previous='Initial', vectors=50)
    mdb.models['Model-1'].steps['Step-1'].setValues(maxIterations=5000)
    
    ## CREATING AND WORKING WITH THE VIRTUAL NODES
    #CREATING SETS
    Part_Full = mdb.models['Model-1'].parts['Part-1']

    Part_Full.Set(name='BOTTOM', nodes=
        Part_Full.nodes.getByBoundingBox(0 - TOL, 0 - TOL, -TOL, numberofunicells*2.*spacing+ TOL, + TOL, TOL) ) 

    Part_Full.Set(name='TOP', nodes=
        Part_Full.nodes.getByBoundingBox(0 - TOL, numberofunicells*2.*spacing - TOL, -TOL, numberofunicells*2.*spacing+ TOL, numberofunicells*2.*spacing+ TOL, TOL) ) 

    Part_Full.Set(name='ALL', edges=
        Part_Full.edges.getByBoundingBox(0 - TOL, 0 - TOL, -TOL, numberofunicells*2.*spacing+ TOL, numberofunicells*2.*spacing+ TOL, TOL) ) 

    Part_Full.Set(name='Left_Nodes', nodes=
        Part_Full.nodes.getByBoundingBox(0 - TOL, 0 - TOL, -TOL, 0 + TOL, numberofunicells*2.*spacing+ TOL, TOL) ) 

    Part_Full.Set(name='Right_Nodes', nodes=
        Part_Full.nodes.getByBoundingBox(numberofunicells*2.*spacing - TOL, 0 - TOL, -TOL, numberofunicells*2.*spacing+ TOL, numberofunicells*2.*spacing+ TOL, TOL) ) 

    # DEFINING MATERIAL PROPERTIES AND SECTION PROPERTIES
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=
        ISOTROPIC, table=((0.5, 0.0), ), testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=POISSON_RATIO)

    #OBTAIN GLOBAL EDGES OF THE ASSEMBLY
    Part_Full = mdb.models['Model-1'].parts['Part-1']
    globalEdgeIndex=[] #all edges
    globalEdgeLIndex=[] #left
    globalEdgeRIndex=[] #righ
    globalEdgeBIndex=[] #bottom
    globalEdgeTIndex=[] #top

    for ii in range(len(Part_Full.edges)):
        pointOnX=Part_Full.edges[ii].pointOn[0][0]
        pointOnY=Part_Full.edges[ii].pointOn[0][1]
        if pointOnX==0.0 or pointOnX==2*spacing*numberofunicells or pointOnY==0.0 or pointOnY==2*spacing*numberofunicells:
            globalEdgeIndex.append(Part_Full.edges[ii:ii+1])
        if pointOnX==0.0:
            globalEdgeLIndex.append(Part_Full.edges[ii:ii+1])
        if pointOnX==2*spacing*numberofunicells:
            globalEdgeRIndex.append(Part_Full.edges[ii:ii+1])
        if pointOnY==0.0:
            globalEdgeBIndex.append(Part_Full.edges[ii:ii+1])
        if pointOnY==2*spacing*numberofunicells:
            globalEdgeTIndex.append(Part_Full.edges[ii:ii+1])


    Part_Full.Set(name='GLOBAL-EDGES', edges=tuple(globalEdgeIndex))

    Part_Full.Set(name='GLOBAL-L-EDGES', edges=tuple(globalEdgeLIndex))
    Part_Full.Set(name='GLOBAL-R-EDGES', edges=tuple(globalEdgeRIndex))
    Part_Full.Set(name='GLOBAL-B-EDGES', edges=tuple(globalEdgeBIndex))
    Part_Full.Set(name='GLOBAL-T-EDGES', edges=tuple(globalEdgeTIndex))

    alteredEdgesIndex=[]
    for ii in range(0,len(Part_Full.edges)):
        pointOnX=Part_Full.edges[ii].pointOn[0][0]
        pointOnY=Part_Full.edges[ii].pointOn[0][1]
        if (pointOnX!=0 and pointOnY!=0 and pointOnY!=2*spacing*numberofunicells and pointOnX!=2*spacing*numberofunicells and ((pointOnX/spacing).is_integer() or (pointOnY/spacing).is_integer())):
            alteredEdgesIndex.append(Part_Full.edges[ii:ii+1])


    Part_Full.Set(name='EDGESALTERED', edges=tuple(alteredEdgesIndex))


    # CREATE REFERENCE POINTS RIGID BODY FOR BOTTOM BASE
    GeometryBottom = (spacing*numberofunicells,0,0)
    GeometryTop = (spacing*numberofunicells,2*spacing*numberofunicells,0)
    
    # CREEATING REFERENCE POINT
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='REF-POINT-BOTTOM', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['REF-POINT-BOTTOM'].ReferencePoint(
        point=GeometryBottom)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='REF-POINT-BOTTOM',
                                                part=mdb.models['Model-1'].parts['REF-POINT-BOTTOM'])
    mdb.models['Model-1'].rootAssembly.Set(name='REF-POINT-BOTTOM', referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],))
    
    # CREATE REFERENCE POINTS RIGID BODY FOR TOP BASE
    # CREEATING REFERENCE POINT
    
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='REF-POINT-TOP', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['REF-POINT-TOP'].ReferencePoint(
        point=GeometryTop)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='REF-POINT-TOP',
                                                part=mdb.models['Model-1'].parts['REF-POINT-TOP'])
    mdb.models['Model-1'].rootAssembly.Set(name='REF-POINT-TOP', referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],))

    # RIGID BODY TIE
    mdb.models['Model-1'].RigidBody(name='BOTTOM-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],
        )), tieRegion=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['BOTTOM'])
    
    mdb.models['Model-1'].RigidBody(name='TOP-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],
        )), tieRegion=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].sets['TOP'])
        
    #APPLY BC
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],)), 
        u1=0.0, u2=0.0, ur3=0.0)
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],)), u1=
        0.0, u2=YDisplacement, ur3=0.0)


    # CREATE PROFILE AND SECTION ASSIGNMENT - THIS IS FOR THE ROUND CROSSECTION
    if DesignB:
        DEDGES = DD
        DDIAGONALS = DEDGES/sqrt(2.)

    elif DesignC:
        DEDGES = DD
        DDIAGONALS = DEDGES/2.

    elif DesignA:
        if Optimization:
            DEDGES=(Lambda/(2.*(1.+Lambda)))*(1.+1./sqrt(2.))/5.
            DDIAGONALS=(1./(2.*2.*(1.+Lambda)))*(1.+sqrt(2.))/5.
            PERIMITERRATIO = PERIMITERRATIO/DEDGES*DD
        else:
            DEDGES = DD
            DDIAGONALS = DEDGES/2.

    else:
        DEDGES = (DD*(1.+sqrt(2.)/2.))
        PERIMITERRATIO = PERIMITERRATIO/(1.+sqrt(2.)/2.)

    #DEFINING SECTION FOR THE GLOBAL EDGES
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES*PERIMITERRATIO, name='GLOBAL-EDGES')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='GLOBAL-EDGES', poissonRatio=0.0, 
        profile='GLOBAL-EDGES', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1'].sets['GLOBAL-EDGES'], sectionName=
        'GLOBAL-EDGES', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1'].sets['GLOBAL-EDGES'])

    # DEFINING SECTION FOR EDGE STRUTS (NON-DIAGONAL)
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES, name='EDGESALTERED')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='EDGESALTERED', poissonRatio=0.0, 
        profile='EDGESALTERED', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1'].sets['EDGESALTERED'], sectionName=
        'EDGESALTERED', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1'].sets['EDGESALTERED'])
    
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
    
    DeleteAbaqusFiles(JobName)
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=THREADS, name=JobName, nodalOutputPrecision=SINGLE, 
        numCpus=8, numDomains=8, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs[JobName].submit(consistencyChecking=OFF)
    mdb.jobs[JobName].waitForCompletion()

    DeleteAbaqusFilesButODB(JobName)
    mdb.saveAs(JobName+'.cae')




def RunPostBuckling(numberofunicells):
    spacing = 1.
    YDisplacement = -0.1*2*spacing*numberofunicells

    Mdb()
    ApplyBuckling(JobName+'.cae',JobName+'.odb',modes,pert,'Step-1','Part-1','Model-1','Part-1-1','ALL')

    mdb.models['Model-1'].StaticStep(initialInc=0.001, maintainAttributes=True, 
        maxNumInc=10000, name='Step-1', nlgeom=ON, previous='Initial')

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-BOTTOM', region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],)), u1=0.0, 
        u2=0.0, ur3=0.0)
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-TOP', region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],)), u1=0.0, 
        u2=YDisplacement, ur3=0.0)

    del mdb.models['Model-1'].boundaryConditions['BC-1']
    del mdb.models['Model-1'].boundaryConditions['BC-2']
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(timeInterval=0.01)

    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=THREADS, name=JobName+'_Post', nodalOutputPrecision=SINGLE, 
        numCpus=8, numDomains=8, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
        userSubroutine='', waitHours=0, waitMinutes=0)

    mdb.jobs[JobName+'_Post'].submit(consistencyChecking=OFF)
    mdb.jobs[JobName+'_Post'].waitForCompletion()
    DeleteAbaqusFilesButODB(JobName+'_Post')



RunBuckling(3,10)
RunPostBuckling(3)

ExtractSetRFDifferentSets('PART-1-1','Step-1',JobName+'_Post','REF-POINT-TOP','REF-POINT-BOTTOM')