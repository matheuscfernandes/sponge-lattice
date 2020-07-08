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

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

def DeleteAbaqusFiles(Job):
	try:
		os.remove(Job+'.odb')
	except: pass
	
	try:
		os.remove(Job+'.dat')
	except: pass
	
	try:
		os.remove(Job+'.com')
	except: pass
	
	try:
		os.remove(Job+'.ipm')
	except: pass
	
	try:
		os.remove(Job+'.log')
	except: pass
	
	try:
		os.remove(Job+'.prt')
	except: pass
	
	try:
		os.remove(Job+'.sim')
	except: pass
	
	try:
		os.remove(Job+'.sta')
	except: pass
		
	try:
		os.remove(Job+'.msg')
	except: pass
	
	try:
		os.remove(Job+'.lck')
	except: pass

	try:
		os.remove(Job+'.inp')
	except: pass

def DeleteAbaqusFilesButODB(Job):
	try:
		os.remove(Job+'.dat')
	except: pass
	
	try:
		os.remove(Job+'.com')
	except: pass
	
	try:
		os.remove(Job+'.ipm')
	except: pass
	
	try:
		os.remove(Job+'.log')
	except: pass
	
	try:
		os.remove(Job+'.prt')
	except: pass
	
	try:
		os.remove(Job+'.sim')
	except: pass
	
	try:
		os.remove(Job+'.sta')
	except: pass
		
	try:
		os.remove(Job+'.msg')
	except: pass
	
	try:
		os.remove(Job+'.lck')
	except: pass
	
	try:
		os.remove(Job+'.inp')
	except: pass

		
		
def ExtractValues(setname,JobName,OutputFile):  
	NameStep='Step-1'
	instancename='PART-1-1'
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[-1]
	Time=odb.steps[NameStep].frames[-1].frameValue
	set=odb.rootAssembly.instances[instancename].nodeSets[setname] 
	displacement=Frame.fieldOutputs['U'].getSubset(region=set).values
	
	for nod in range(0,len(set.nodes)):
				xDisp=displacement[nod].data[0]
				yDisp=displacement[nod].data[1]
				zDisp=displacement[nod].data[2]
				
				xCoor=set.nodes[nod].coordinates[0]
				yCoor=set.nodes[nod].coordinates[1]
				
				
				text_file.write('%e %e %e %e %e %e\r\n' % (Time, xCoor, yCoor, xDisp, yDisp, zDisp))
	text_file.close()
			
	
def ReturnRF(instancename,setname,JobName):  
	NameStep='Step-1'
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[FrameLength-1]
	Time=odb.steps[NameStep].frames[FrameLength-1].frameValue
	set=odb.rootAssembly.nodeSets[setname] 
	displacement=Frame.fieldOutputs['RF'].getSubset(region=set).values
	xRF=[]
	yRF=[]
	zRF=[]

	xCoor=[]
	yCoor=[]
	for nod in range(0,len(set.nodes)):
		xRF.append(displacement[nod].data[0])
		yRF.append(displacement[nod].data[1])
		zRF.append(displacement[nod].data[2])

	odb.close()
	return np.sum(xRF), np.sum(yRF), np.sum(zRF)

def ExtractRF2D(setname,JobName,OutputFile):  
	NameStep='Step-1'
	instancename='PART-1-1'
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[NameStep].frames)
	Frame=odb.steps[NameStep].frames[FrameLength-1]
	Time=odb.steps[NameStep].frames[FrameLength-1].frameValue
	set=odb.rootAssembly.instances[instancename].nodeSets[setname] 
	displacement=Frame.fieldOutputs['RF'].getSubset(region=set).values
	
	for nod in range(0,len(set.nodes)):
				xRF=displacement[nod].data[0]
				yRF=displacement[nod].data[1]
				
				xCoor=set.nodes[nod].coordinates[0]
				yCoor=set.nodes[nod].coordinates[1]
				
				
				text_file.write('%e %e %e %e %e\r\n' % (Time, xCoor, yCoor, xRF, yRF))	
	text_file.close()
	odb.close()

def ExtractVirtualPointRF(instanceVP,nodeSetVP,stepName,JobName,OutputFile):
	text_file = open(OutputFile, 'w')
	odb=openOdb(path=JobName+'.odb')
	FrameLength=len(odb.steps[stepName].frames)
	set = odb.rootAssembly.instances[instanceVP].nodeSets[nodeSetVP]
	for i in range(FrameLength):
		Frame = odb.steps[stepName].frames[i].frameValue
		RF = odb.steps[stepName].frames[i].fieldOutputs['RF'].getSubset(region = set).values[0].data[1]
		text_file.write('%e %e\r\n' % (Frame, RF))
	text_file.close()
	odb.close()

def ApplyBuckling(mdbName,odbName,ImpFrames,ImpWeights,StepName,PartName,ModelName,InstanceName,AllNodes):
	########## DESCRIPTION OF VARIABLES FOR THE DEFINED FUNCTION
	'''
	mdbName -- 	THIS IS THE NAME FOR THE MDB PATH FOR THE CAE
	odbName -- THIS IS THE PATH FOR THE ODB FILE
	ImpFrames -- THIS IS THE FRAME NUMBERS OF THE FREQUENCY NUMBER OF THE ANALYSIS
	ImpWeigts -- THIS IS THE WEIGHT FOR THE FREQUENCY DISPLACEMENT B.C.
	StepName -- THIS IS THE NAME OF THE STEP
	PartName -- THIS IS THE NAME OF THE PART WHERE THE BUCKLING WILL BE APPLIED
	ModelName -- THIS IS THE MODEL NAME IN THE CAE
	InstanceName -- THIS IS THE INSTANCE NAME IN THE CAE
	AllNodes -- THIS IS THE SET NAME THAT CONTAINS ALL OF THE NODES
	'''
	##########
	# BEGIN CODE BY IMPORTING THE MDB MODEL OF FREQUENCY ANALYSIS
	mdb=openMdb(pathName=mdbName)
	odb=openOdb(path=odbName)
	# THIS OPENS THE RESULTS ODB WITH THE NODESET OF ALL OF THE NODES
	pbpPartNodes=odb.rootAssembly.instances[InstanceName.upper()].nodeSets[AllNodes.upper()]

	#CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
	NewCoord=np.zeros((len(mdb.models[ModelName].parts[PartName].nodes), 3))	
	
	for CImp in range(len(ImpFrames)):
		cframe = ImpFrames[CImp]
		firstFrame = odb.steps[StepName].frames[cframe]
		displacement = firstFrame.fieldOutputs['U']
		pbpDispField = displacement.getSubset(region=pbpPartNodes)
		pbpDisp = pbpDispField.values


		# Imperfection Using Buckling Analysis Results
		#---------------------------------------------------------------
		ind=0;
		IMP = ImpWeights[CImp]
		#CALCULATE THE MODIFIED COORDINATES
		for i in mdb.models[ModelName].parts[PartName].nodes:
			NewCoord[ind][0]=i.coordinates[0]+IMP*pbpDisp[ind].data[0]
			NewCoord[ind][1]=i.coordinates[1]+IMP*pbpDisp[ind].data[1]
			#NewCoord[ind][2]=i.coordinates[2]+IMP*pbpDisp[ind].data[2]
			ind=ind+1

		#SET THE NEW COORDINATES
		mdb.models[ModelName].parts[PartName].editNode(
			nodes=mdb.models[ModelName].parts[PartName].nodes,
			coordinates=NewCoord)

	mdb.models[ModelName].rootAssembly.regenerate()

	print 'Original Coordinates Modified Successfully!'

def ExtractSetDisplacement(instance,stepName,JobName,nodeSet1):
	text_file = open(JobName+'_DispOutput.txt', 'w')
	odb=openOdb(path=JobName+'.odb')
	set1 = odb.rootAssembly.instances[instance].nodeSets[nodeSet1]
	allNodes=[]
	for node in odb.rootAssembly.instances[instance].nodes:
		allNodes.append([node.label,node.coordinates])

	allNodes=sorted(allNodes, key=lambda x : x[0])

	frameNum=0
	for frame in odb.steps[stepName].frames:
		i=0
		for s1 in frame.fieldOutputs['U'].getSubset(region = set1).values:

			xdisp=s1.data[0]
			ydisp=s1.data[1]

			coor=allNodes[s1.nodeLabel-1][1]
			text_file.write('{} {} {} {} {} {}\r\n'.format(coor[0],coor[1],coor[2],frameNum,xdisp,ydisp))
			i+=1
		frameNum+=1


	text_file.close()
	
	odb.close()

def ExtractSetRF(instance,stepName,JobName,nodeSet1):
	text_file = open(JobName+'_Output.txt', 'w')
	odb=openOdb(path=JobName+'.odb')
	set1 = odb.rootAssembly.instances[instance].nodeSets[nodeSet1]
	allNodes=[]
	for node in odb.rootAssembly.instances[instance].nodes:
		allNodes.append([node.label,node.coordinates])

	allNodes=sorted(allNodes, key=lambda x : x[0])

	for frame in odb.steps[stepName].frames:
		for disp,rf in zip(frame.fieldOutputs['U'].getSubset(region = set1).values,frame.fieldOutputs['RF'].getSubset(region = set1).values):

			xdisp=disp.data[0]
			ydisp=disp.data[1]

			xrf=rf.data[0]
			yrf=rf.data[1]

			coor=allNodes[disp.nodeLabel-1][1]
			text_file.write('{} {} {} {} {} {} {}\r\n'.format(coor[0],coor[1],coor[2],xdisp,ydisp,xrf,yrf))

	text_file.close()
	
	odb.close()
