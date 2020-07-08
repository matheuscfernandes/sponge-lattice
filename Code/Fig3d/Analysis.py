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
import meshEdit


session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

execfile('Functions.py')


def RunSimulation(numberofunicells,structurenumberofunicells,THICKNESS,Angle):
    Mdb()

    PERIMITERRATIO=2.
    spacing = 1
    TOL = 10E-6
    DD = 0.1*spacing*1.1
    AppYLoad=-1.
    x2 = spacing / (sqrt(2) + 2)

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
        #x1=sqrt(2)*spacing/(sqrt(2)*2)
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

    HSize=structurenumberofunicells*spacing
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)

    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-HSize, -HSize), point2=(-HSize,HSize))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(-HSize, HSize), point2=(HSize,HSize))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(HSize, HSize), point2=(HSize,-HSize))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(HSize, -HSize), point2=(-HSize,-HSize))


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

    #CREATING INSTANCE AND MESHING
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    Instant_Full=mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])

    #CREATING LINEAR PATTERN FOR FINITE SIZE
    mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
        0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Part-1-1', ), number1=numberofunicells, 
        number2=numberofunicells, spacing1=2.*spacing, spacing2=2.*spacing)
    
    #move the part
    mdb.models['Model-1'].rootAssembly.translate(instanceList=mdb.models['Model-1'].rootAssembly.allInstances.keys(), 
        vector=(-numberofunicells*spacing, -numberofunicells*spacing, 0.0))

    mdb.models['Model-1'].rootAssembly.rotate(angle=Angle, axisDirection=(0.0, 0.0, 
        1.0), axisPoint=(0.0, 0.0, 0.0), instanceList=mdb.models['Model-1'].rootAssembly.allInstances.keys())

    #ADD OUTER EDGE PART 2 INTO THE ASSEMBLY
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
        part=mdb.models['Model-1'].parts['Part-2'])
    
    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=BOTH, 
        instances=(mdb.models['Model-1'].rootAssembly.instances.values()), mergeNodes=
        BOUNDARY_ONLY, name='Part-3', nodeMergingTolerance=1e-09, 
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

    Part_Full = mdb.models['Model-1'].parts['Part-1']
    Part_Full.Set(name='OfInterest', elements=
        Part_Full.elements.getByBoundingBox(-structurenumberofunicells*spacing - TOL, -structurenumberofunicells*spacing - TOL, -TOL, structurenumberofunicells*spacing+ TOL, structurenumberofunicells*spacing+ TOL, TOL) ) 

    mdb.models['Model-1'].parts['Part-1'].PartFromMesh(copySets=True, name='Part-1-mesh-1')

    mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
    Part_Full = mdb.models['Model-1'].parts['Part-1-mesh-1']
    
    edgesPlanes=structurenumberofunicells*spacing    
    maxPlane=2*spacing*numberofunicells

    def DelteElements(X1,Y1,X2,Y2):
        delElements=Part_Full.elements.getByBoundingBox(X1 + TOL, Y1 + TOL, -TOL, X2 - TOL, Y2 - TOL, TOL)
        Part_Full.deleteElement(elements=delElements)

    DelteElements(-maxPlane,-maxPlane,-edgesPlanes,maxPlane)
    DelteElements(-maxPlane,-maxPlane,maxPlane,-edgesPlanes)
    DelteElements(edgesPlanes,-maxPlane,maxPlane,maxPlane)
    DelteElements(-maxPlane,edgesPlanes,maxPlane,maxPlane)

    #CREATING STEP
    mdb.models['Model-1'].BuckleStep(name='Buckle',
                                   numEigen=6, previous='Initial', vectors=50)
    mdb.models['Model-1'].steps['Buckle'].setValues(maxIterations=5000)

    ## CREATING AND WORKING WITH THE VIRTUAL NODES
    #CREATING SETS
    Part_Full = mdb.models['Model-1'].parts['Part-1-mesh-1']

    HSize=structurenumberofunicells*spacing


    Part_Full.Set(name='BOTTOM', nodes=
        Part_Full.nodes.getByBoundingBox(-HSize - TOL, -HSize - TOL, -TOL,HSize+ TOL, -HSize+ TOL, TOL) ) 

    Part_Full.Set(name='TOP', nodes=
        Part_Full.nodes.getByBoundingBox(-HSize - TOL, HSize - TOL, -TOL,HSize+ TOL, HSize+ TOL, TOL) )
        
    Part_Full.Set(name='LEFT', nodes=
        Part_Full.nodes.getByBoundingBox(-HSize - TOL, -HSize - TOL, -TOL, -HSize+ TOL, HSize+ TOL, TOL) )

    Part_Full.Set(name='RIGHT', nodes=
        Part_Full.nodes.getByBoundingBox(HSize - TOL, -HSize - TOL, -TOL,HSize+ TOL, HSize+ TOL, TOL) )



    Part_Full.Set(name='EDGES-PART2', elements=
        [Part_Full.elements.getByBoundingBox(-HSize - TOL, -HSize - TOL, -TOL,HSize+ TOL, -HSize+ TOL, TOL) ,
        Part_Full.elements.getByBoundingBox(-HSize - TOL, HSize - TOL, -TOL,HSize+ TOL, HSize+ TOL, TOL),
        Part_Full.elements.getByBoundingBox(-HSize - TOL, -HSize - TOL, -TOL, -HSize+ TOL, HSize+ TOL, TOL),
        Part_Full.elements.getByBoundingBox(HSize - TOL, -HSize - TOL, -TOL,HSize+ TOL, HSize+ TOL, TOL)])

    mdb.models['Model-1'].rootAssembly.regenerate()

    
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=
        'Part-1-mesh-1-1', part=mdb.models['Model-1'].parts['Part-1-mesh-1'])
    del mdb.models['Model-1'].rootAssembly.features['Part-1-1']
    
    # DEFINING MATERIAL PROPERTIES AND SECTION PROPERTIES
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Hyperelastic(materialType=
        ISOTROPIC, table=((0.5, 0.0), ), testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=POISSON_RATIO)
    

    # CREATE PROFILE AND SECTION ASSIGNMENT - THIS IS FOR THE ROUND CROSSECTION
    if DesignB:
        DEDGES = DD
        DDIAGONALS = DEDGES/sqrt(2.)

    elif DesignC:
        DEDGES = DD
        DDIAGONALS = DEDGES/2.

    elif DesignA:
        DEDGES = DD
        DDIAGONALS = DEDGES/2.

    else:
        DEDGES = (DD*(1.+sqrt(2.)/2.))
        PERIMITERRATIO = PERIMITERRATIO/(1.+sqrt(2.)/2.)


    # DEFINING SECTION FOR EDGE STRUTS (NON-DIAGONAL)
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES, name='EDGES')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='EDGES', poissonRatio=0.0, 
        profile='EDGES', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['EDGES'], sectionName=
        'EDGES', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['EDGES'])

    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES, name='MIDDLE')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='MIDDLE', poissonRatio=0.0, 
        profile='MIDDLE', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['MIDDLE'], sectionName=
        'MIDDLE', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['MIDDLE'])
    
    # DEFINING SECTION FOR THE DIAGONAL STRUTS
    if DesignB or DesignA or DesignC:
        mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DDIAGONALS, name='DIAGONAL')
        mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
            DURING_ANALYSIS, material='Material-1', name='DIAGONAL', poissonRatio=0.0, 
            profile='DIAGONAL', temperatureVar=LINEAR)
        mdb.models['Model-1'].parts['Part-1-mesh-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=
            mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['DIAGONAL'], sectionName=
            'DIAGONAL', thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['Part-1-mesh-1'].assignBeamSectionOrientation(method=
            N1_COSINES, n1=(0.0, 0.0, -1.0), region=
            mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['DIAGONAL'])

    #DEFINING SECTION FOR THE GLOBAL EDGES
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=DEDGES*PERIMITERRATIO, name='EDGES-PART2')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=
        DURING_ANALYSIS, material='Material-1', name='EDGES-PART2', poissonRatio=0.0, 
        profile='EDGES-PART2', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['EDGES-PART2'], sectionName=
        'EDGES-PART2', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1-mesh-1'].assignBeamSectionOrientation(method=
        N1_COSINES, n1=(0.0, 0.0, -1.0), region=
        mdb.models['Model-1'].parts['Part-1-mesh-1'].sets['EDGES-PART2'])
    
    # CREATE REFERENCE POINTS RIGID BODY FOR BOTTOM BASE
    GeometryBottom = (0,-structurenumberofunicells*spacing,0)
    GeometryTop = (0,structurenumberofunicells*spacing,0)
    
    # CREEATING REFERENCE POINT
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='REF-POINT-BOTTOM', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['REF-POINT-BOTTOM'].ReferencePoint(
        point=GeometryBottom)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='REF-POINT-BOTTOM',
                                                part=mdb.models['Model-1'].parts['REF-POINT-BOTTOM'])
    mdb.models['Model-1'].rootAssembly.Set(name='REF-POINT-BOTTOM', referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],))
    
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
    Part_Full = mdb.models['Model-1'].rootAssembly.instances['Part-1-mesh-1-1']
    mdb.models['Model-1'].RigidBody(name='BOTTOM-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],
        )), tieRegion=Part_Full.sets['BOTTOM'])
    
    mdb.models['Model-1'].RigidBody(name='TOP-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],
        )), tieRegion=Part_Full.sets['TOP'])
    
    #APPLY BC
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, createStepName='Buckle', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=
        Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],)), 
        u1=0.0, u2=0.0, ur3=0.0)

    mdb.models['Model-1'].ConcentratedForce(cf1=0.0, cf2=AppYLoad, createStepName=
            'Buckle', distributionType=UNIFORM, field='', localCsys=None,  name='BC-REF-2', 
            region=Region(referencePoints=(mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],)))

    # CREATE JOB AND RUN IT
    if DesignB:
        JobName="DesignB"
    elif DesignA:
        JobName="DesignA"
    elif DesignC:
        JobName="DesignC"
    else:
        JobName="DesignD"
    
    FileName=JobName+'_Output.txt'

    JobName=JobName+'_'+str(int(Angle))
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
    
    EV = np.array(ExtractEigenMode(JobName, 6))*AppYLoad
    FileWrite=open(FileName, 'a')
    FileWrite.write('{} {} {} {} {} {} {}\r\n'.format(Angle,EV[0],EV[1],EV[2],EV[3],EV[4],EV[5]))
    FileWrite.close()
    if Angle not in [0,10,20,30,40,45]:
        DeleteAbaqusFiles(JobName)

# allAngles=np.linspace(0,1,91)*90

for Angle in allAngles:

    RunSimulation(15,10,10,Angle)

