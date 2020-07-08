execfile('Analysis.py')
execfile('Parameters.py')
execfile('Global_Parameters.py')


Return=RunModel(Lambda,totalDiags,AllDiagsSeparation)

wfile = open('Return_'+str(WorkerNumber)+'.txt', 'w')
wfile.write('%e' % (Return))
wfile.close()