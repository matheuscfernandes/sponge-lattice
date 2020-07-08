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

from Functions import *

import numpy as np
import math
import time	

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

####### MODEL FUNCTIONS #######

TOL = 10E-6  # periodic boundary search tolerance
AppYLoad = -1.  # applied y load boundary condition for buckling portion

YoungsModulus = 1.0 # this is set to scale the stress values to become closer to values seen critical buckling

PoissonsRatio = 0.3
THICKNESS = 100.

numRep=3

def RunModel(Lambda,totalDiags,AllDiagsSpacing):
    Mdb()

    diagNum=totalDiags

    # LOGISTICS FOR DECIDING THE NUMBER OF DIAGONALS DURING GEOMETRY CONSTRUCTION
    if totalDiags==1:
        oddNumberDiags=True
        evenNumberDiags=False

    elif totalDiags%2==1:
        totalDiags=(totalDiags-1)
        oddNumberDiags=True
        evenNumberDiags=True
    else:
        oddNumberDiags=False
        evenNumberDiags=True

    # CREATING THE MAIN UNICEL
    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__', sheetSize=200.0)

    # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, spacing / 2.), point2=(
    #     2. * spacing, spacing / 2.))
    # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 3. * spacing / 2.), point2=(
    #     2. * spacing, 3. * spacing / 2.))

    # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing / 2., 0), point2=(
    #     spacing / 2., 2. * spacing))
    # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(3. * spacing / 2., 0), point2=(
    #     3. * spacing / 2., 2. * spacing))


    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
        2. * spacing, 0.0))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, spacing), point2=(
        2. * spacing, spacing ))

    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 0), point2=(
        0, 2. * spacing))
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
        spacing, 2. * spacing))

    # CONSTRUCTION OF ODD CENTER DIAGONAL REINFORCEMENT
    if oddNumberDiags:
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.0, spacing), point2=(
            spacing, 2 * spacing))
        # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 2 * spacing), point2=(
        #     2 * spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0 , 2 * spacing), point2=(
            spacing, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
            2 * spacing, spacing))
        # mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spacing, 0), point2=(
        #     0, spacing))
        mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing, 0), point2=(
            spacing, spacing))

    # CONSTRUCTION OF THE EVEN NUMBER OF DIAGONAL REINFORCEMENT
    if evenNumberDiags:
        for DiagNum in range(totalDiags/2):
            # CREATING UNICEL FOR THE TWO DIAGONAL PART
            DiagSeparation=AllDiagsSpacing[DiagNum]

            if DiagSeparation<0.001:
                print('Diag Separation too small -- corrected to 0 diag separation value')
                DiagSeparation=0
            if DiagSeparation>spacing*0.999:
                print('Diag Separation too large -- corrected to maximum values')
                DiagSeparation=spacing
                
            shift =  DiagSeparation

            sp1 = shift
            sp2 = spacing - shift
            sp3 = spacing + shift
            sp4 = 2* spacing -shift

            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, sp1), point2=(
                sp1, 0))
            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(2*spacing, sp4), point2=(
                sp4, 2*spacing))

            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, sp2), point2=(
                sp3,2*spacing))
            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, sp3), point2=(
                sp2,2*spacing))

            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, sp4), point2=(
                sp4,0))
            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(sp1,2*spacing), point2=(
                2*spacing,sp1))

            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(sp2, 0), point2=(
                2*spacing,sp3))
            mdb.models['Model-1'].sketches['__profile__'].Line(point1=(sp3,0), point2=(
                2*spacing,sp2))

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                            name='Part-1', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseWire(
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    # CREATING OUTER EDGE PART
    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0, 0), point2=(
        2*spacing*numRep, 0))
    
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                            name='Part-Edge', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-Edge'].BaseWire(
        sketch=mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']

    # CREATING STEP
    mdb.models['Model-1'].BuckleStep(name='Step-1',
                                 numEigen=8, previous='Initial', vectors=28, maxIterations=10000)

    # SECTION DEFINITIONS FOR THE DIFFERENT PARTS OF THE STRUCTURE
    mdb.models['Model-1'].rootAssembly.regenerate()
    Part_Full = mdb.models['Model-1'].parts['Part-1']
    edgeIndex = []
    diagIndex = []
    for ii in range(len(Part_Full.edges)):
        pointOnX = Part_Full.edges[ii].pointOn[0][0]
        pointOnY = Part_Full.edges[ii].pointOn[0][1]
        if pointOnX == 0 or pointOnX == spacing or pointOnY == 0 or pointOnY == spacing:
            edgeIndex.append(Part_Full.edges[ii:ii + 1])
        else:
            diagIndex.append(Part_Full.edges[ii:ii + 1])
    edgeIndex = tuple(edgeIndex)
    diagIndex = tuple(diagIndex)
    Part_Full.Set(name='DIAGONAL', edges=diagIndex)
    Part_Full.Set(name='EDGES', edges=edgeIndex)


    Part_Full = mdb.models['Model-1'].parts['Part-Edge']
    Part_Full.Set(name='OUT-EDGES', edges=Part_Full.edges[:])

    # DEFINING MATERIAL PROPERTIES AND SECTION PROPERTIES
    mdb.models['Model-1'].Material(name='Material-1')
    mdb.models['Model-1'].materials['Material-1'].Elastic(
        table=((YoungsModulus, PoissonsRatio), ))

    # COMPUTING THICKNESSES OF MEMBERS BASED ON VOLUME RATIO, LAMBDA AND MASS CONSERVATION
    # TND = D*(1+sqrt(2)/2)*(1+1/Lambda)**-1
    # TD  = D/sqrt(2)*(1+sqrt(2))*(1+Lambda)**-1

    TND=(Lambda/(2.*(1.+Lambda)))*(1.+1./sqrt(2.))
    TD=(1./(2.*diagNum*(1.+Lambda)))*(1.+sqrt(2.))

    # DEFINING SECTION FOR EDGE STRUTS (NON-DIAGONAL)
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=TND, name='EDGES')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS, material='Material-1', name='EDGES', poissonRatio=0.0,
                                    profile='EDGES', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0,
                                                            offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models['Model-1'].parts['Part-1'].sets['EDGES'], sectionName='EDGES', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=N1_COSINES, n1=(
        0.0, 0.0, -1.0), region=mdb.models['Model-1'].parts['Part-1'].sets['EDGES'])

    # DEFINING SECTION FOR THE DIAGONAL STRUTS
    mdb.models['Model-1'].RectangularProfile(a=THICKNESS, b=TD, name='DIAGONAL')
    mdb.models['Model-1'].BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS, material='Material-1', name='DIAGONAL', poissonRatio=0.0,
                                    profile='DIAGONAL', temperatureVar=LINEAR)
    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0,
                                                            offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models['Model-1'].parts['Part-1'].sets['DIAGONAL'], sectionName='DIAGONAL', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=N1_COSINES, n1=(
        0.0, 0.0, -1.0), region=mdb.models['Model-1'].parts['Part-1'].sets['DIAGONAL'])

    # DEFINING SECTION FOR OUTER EDGE STRUTS (NON-DIAGONAL)
    mdb.models['Model-1'].parts['Part-Edge'].SectionAssignment(offset=0.0,
                                                            offsetField='', offsetType=MIDDLE_SURFACE, region=mdb.models['Model-1'].parts['Part-Edge'].sets['OUT-EDGES'], sectionName='EDGES', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].parts['Part-1'].assignBeamSectionOrientation(method=N1_COSINES, n1=(
        0.0, 0.0, -1.0), region=mdb.models['Model-1'].parts['Part-Edge'].sets['OUT-EDGES'])    

    # CREATING INSTANCE
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1',
        part=mdb.models['Model-1'].parts['Part-1'])
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.LinearInstancePattern(direction1=(1.0, 0.0, 
        0.0), direction2=(0.0, 1.0, 0.0), instanceList=('Part-1-1', ), number1=numRep, 
        number2=numRep, spacing1=2*spacing, spacing2=2*spacing)

    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-Edge-1',
        part=mdb.models['Model-1'].parts['Part-Edge'])
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-Edge-2',
        part=mdb.models['Model-1'].parts['Part-Edge'])
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-Edge-2', ), 
        vector=(0.0, 2*spacing*numRep, 0.0))
    mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 0.0, 
        1.0), axisPoint=(2*spacing*numRep, 0.0, 0.0), instanceList=('Part-Edge-1', ))

    mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
        instances=(mdb.models['Model-1'].rootAssembly.instances.values()), name=
        'Part-2', originalInstances=SUPPRESS)
    mdb.models['Model-1'].rootAssembly.regenerate()

    # MESHING FINAL PART
    mdb.models['Model-1'].parts['Part-2'].seedPart(deviationFactor=0.25,
    minSizeFactor=0.25, size=0.25)

    mdb.models['Model-1'].parts['Part-2'].Set(
        name='ALL', edges=mdb.models['Model-1'].parts['Part-2'].edges[:])
    mdb.models['Model-1'].parts['Part-2'].setElementType(elemTypes=(ElemType(
        elemCode=B22, elemLibrary=STANDARD), ), regions=mdb.models['Model-1'].parts['Part-2'].sets['ALL'])
    mdb.models['Model-1'].parts['Part-2'].generateMesh()
    mdb.models['Model-1'].rootAssembly.regenerate()

    # CREATING SETS
    Instance_Full = mdb.models['Model-1'].rootAssembly.instances['Part-2-1']

    mdb.models['Model-1'].rootAssembly.Set(name='TOP', nodes=Instance_Full.nodes.getByBoundingBox(
        0.-TOL, 2. * spacing*numRep - TOL, 0-TOL, 2. * spacing*numRep + TOL, 2. * spacing*numRep + TOL, 0.+TOL))

    mdb.models['Model-1'].rootAssembly.Set(name='BOTTOM', nodes=Instance_Full.nodes.getByBoundingBox(
        0.-TOL, 0 - TOL, 0-TOL, 2. * spacing*numRep + TOL, 0 + TOL, 0.+TOL))

    # CREEATING REFERENCE POINT FOR BOTTOM
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='REF-POINT-BOTTOM', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['REF-POINT-BOTTOM'].ReferencePoint(
        point=(spacing*numRep,0,0))
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='REF-POINT-BOTTOM',
                                                part=mdb.models['Model-1'].parts['REF-POINT-BOTTOM'])
    mdb.models['Model-1'].rootAssembly.Set(name='REF-POINT-BOTTOM', referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],))
    
    # CREEATING REFERENCE POINT FOR TOP

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='REF-POINT-TOP', type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['REF-POINT-TOP'].ReferencePoint(
        point=(spacing*numRep,2*spacing*numRep,0))
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='REF-POINT-TOP',
                                                part=mdb.models['Model-1'].parts['REF-POINT-TOP'])
    mdb.models['Model-1'].rootAssembly.Set(name='REF-POINT-TOP', referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],))

    # RIGID BODY TIE
    mdb.models['Model-1'].RigidBody(name='BOTTOM-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1],
        )), tieRegion=mdb.models['Model-1'].rootAssembly.sets['BOTTOM'])
    
    mdb.models['Model-1'].RigidBody(name='TOP-CONSTRAINT', refPointRegion=Region(
        referencePoints=(
            mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1],
        )), tieRegion=mdb.models['Model-1'].rootAssembly.sets['TOP'])


    # APPLY BOUNDARY CONDITIONS
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-BOTTOM-STEP-0', region=Region(referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-BOTTOM'].referencePoints[1], 
        )), u1=0.0, u2=0.0, ur3=0.0)

    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-TOP-STEP-0', region=Region(referencePoints=(
        mdb.models['Model-1'].rootAssembly.instances['REF-POINT-TOP'].referencePoints[1], 
        )), u1=0.0, u2=UNSET, ur3=0.0)

    mdb.models['Model-1'].ConcentratedForce(cf2=AppYLoad, createStepName='Step-1', 
        distributionType=UNIFORM, field='', localCsys=None, name='Load', region=
        mdb.models['Model-1'].rootAssembly.sets['REF-POINT-TOP'])

    # DeleteAbaqusFiles(JobName)
    
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
            memory=6, memoryUnits=GIGA_BYTES, model='Model-1', modelPrint=OFF,
            multiprocessingMode=DEFAULT, name=JobName, nodalOutputPrecision=SINGLE,
            numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS,
            userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs[JobName].submit(consistencyChecking=OFF)
    mdb.jobs[JobName].waitForCompletion()
    
    # EXTRACT INFORMATION INTO FILES
    EigenValues = np.array(ExtractEigenMode(JobName, 8))*(-1.*AppYLoad) 

    # DeleteAbaqusFiles(JobName)

    return (EigenValues[0]) 
