import subprocess
import os
from scipy.optimize import minimize
import numpy as np
import cma
import random
import shutil

exec(open("./Global_Parameters.py").read())

#Clean created directory folder by changing directory and deleting the folder
def CleanWorkerDirectory(WorkerFolder,owd):
	os.chdir(owd)
	shutil.rmtree(WorkerFolder)

# Create directory folder and change directory
def CreateWorkerDirectory(WorkerNumber):
	#name of folder
	WorkerFolder='Files_'+str(WorkerNumber)

	#obtain current path
	owd = os.getcwd()

	#create new folder with above defined name
	os.mkdir(WorkerFolder)
	
	#copy common files for simulations into new directory
	shutil.copy2('Functions.py',WorkerFolder)
	shutil.copy2('Functions.pyc',WorkerFolder)
	shutil.copy2('Analysis.py',WorkerFolder)
	shutil.copy2('Global_Parameters.py',WorkerFolder)

	#move specific files for simulatiosn into new directory
	shutil.move('MainRun_'+str(WorkerNumber)+'.py',WorkerFolder)
	shutil.move('Parameters_'+str(WorkerNumber)+'.py',WorkerFolder)

	#change current directory into new directory
	os.chdir(WorkerFolder)

	#return the new directory folder name and the old directory path
	return WorkerFolder, owd


#Optimization objective function
def objective(X):
	#Extracting parameters for passing to simulations
	Lambda=X[0]*1000 #normalizing by 1000
	AllDiagsSeparation=list(np.array([ X[1],X[2],X[3],X[4] ])*spacing) #normalizing by spacing

	#creating a random worker number
	WorkerNumber=random.randint(1,9999999999)
	
	#creating a spefific files with specific parameters for a given worker
	wfile = open('Parameters_'+str(WorkerNumber)+'.py', 'w')
	wfile.write('Lambda = {}\r\n'.format(Lambda))
	wfile.write('AllDiagsSeparation = '+str(AllDiagsSeparation)+'\r\n')
	wfile.write('WorkerNumber = '+str(WorkerNumber)+'\r\n')
	wfile.write('JobName = "DesignA_DiagNum"+str('+str(totalDiags)+')+"_"+str(WorkerNumber)\r\n')
	wfile.write('FileName = JobName+"_Output.txt"'+'\r\n')
	wfile.close()

	#open the common main run file to edit to open specific parameters file
	with open('MainRun.py', 'r') as file:
		# read a list of lines into data
		MainRunFile = file.readlines()

	# now change the 2nd line, note that you have to add a newline
	MainRunFile[1] = 'execfile(\'Parameters_'+str(WorkerNumber)+'.py\')\r\n'

	# and write everything back
	with open('MainRun_'+str(WorkerNumber)+'.py', 'w') as file:
		file.writelines( MainRunFile )

	# create and move files and change directory to specific worker directory
	# this step is important for abaqus to avoid using similar meta-files at the same time - without this step abaqus breaks 
	WorkerFolder, owd=CreateWorkerDirectory(WorkerNumber)
	
	# run abaqus cae using a python script
	FNULL = open(os.devnull, 'w')
	subprocess.call('abaqus cae noGUI=MainRun_'+str(WorkerNumber)+'.py',shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
	FNULL.close()

	# obtain final value from the output of abaqus
	rfile = open('Return_'+str(WorkerNumber)+'.txt', 'r')
	returnval = float(rfile.read())
	rfile.close()
	
	# remove worker directory and all specific worker files
	CleanWorkerDirectory(WorkerFolder,owd)

	#return the inverse of the critical buckling strain -- since we want to maximize strain we want to minimize the inverse of strain
	return 1/returnval

#main optimization file. It is important to assure the name is main such that each worker also does not run this section in parallel
if __name__=='__main__':
	
	X0=[1.]+[0.5]*4
	sigma0=0.1

	# Initialize the CMA-ES optimizer (All variables lie in the interval 0,1)
	es = cma.CMAEvolutionStrategy(X0, sigma0, {'bounds': [0.0, 1.0], 'maxiter': 200,'popsize':20})
		
	# EvalParallel() enables the paralle evaluation of candidate points every generation
	# This greatly speeds up the computation time for each optimization run
	with cma.fitness_transformations.EvalParallel(20) as parallel:
		
		# While not optimal keep iterating
		while not es.stop():
			
			# Ask for a new solution
			X = es.ask()

			# Evaluate the solution using the parallel solver
			es.tell(X, parallel(objective, X))

			# Display progress on screen
			es.disp(1)

	wfile = open('FinalParamsAndVal_DiagNum'+str(totalDiags)+'.txt', 'w')
	wfile.write(str(es.result))
	wfile.close()
	# Result documentation: http://cma.gforge.inria.fr/apidocs-pycma/cma.evolution_strategy.CMAEvolutionStrategyResult.html