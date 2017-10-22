#########################
### INITIALISE PYTHON ###
#########################

#execfile('in_plots.py')					# initialise python to allow plotting
#execfile('in.py')						# initialise python

#set_printoptions(threshold='nan')
#set_printoptions(linewidth=1000)

import os
import numpy
from numpy import *

import scipy
from scipy import *

#import matplotlib
#matplotlib.use('macosx')
#from matplotlib import *
#import matplotlib.pyplot as plt

#import pylab
#from pylab import *

import random
import time


###########################
### CONSTANT PARAMETERS ###
###########################

# world parameters
x_min = -4000
x_max = 4000
y_min = -4000
y_max = 4000

# population parameters
generations = 20001											# number of steps in random walk
in_entities = 1												# initial number of entities
max_entities = 10000										# carrying capacity - maximum number of entities the modelled domain can sustain

# craniometric (CM) measurement parameters
num_CM_meas = 50											# total number of CM measurements considered

# gaussian distribution parameters - for CM measurements
CM_mu = 1000 												# mean
CM_sigma_frac = 0.05											# standard deviation of CM gaussian - proportion of mean (CM_mu): 0.005, 0.05, 0.5
CM_sigma_val = CM_mu*CM_sigma_frac							# standard deviation of CM gaussian - value (CM_mu*CM_sigma_frac)

prob_fis = 0.1 											# probability with which each entity fissions at each generation
prob_dep = 0.00001											# probability with which each entity is sampled (data recorded to output)

for taskID in range(1,2):
	
	###########################
	### VARIABLE PARAMETERS ###
	###########################
	
	# gaussian distribution parameters - for migration
	mig_mu = 0													# mean
#	mig_sigma = random.randint(1,125)	# 2.5 #					# standard deviation: 2.5, 5, 7.5, 10
	mig_sigma = random.uniform(0,100)
	
	##################################################
	### SAVE MODEL PARAMETER VALUES TO OUTPUT FILE ###
	##################################################
	
	### record parameters to file
	##q = open('sim_data_2D_mig{0}_CM{1}.txt'.format("%.2f"% round(mig_sigma/float(25),2),CM_sigma_frac),'a')
	##q.write("PARAMS VARIED BETWEEN SIMS:"+"\n"+'average distance (km/yr):     '+str(mig_sigma/float(25))+"\n"+"\n"+"PARAMS FIXED BETWEEN SIMS:"+'\n'+"generations:                  "+str(generations)+'\n'+'\n'+'number of CM measurements:    '+str(num_CM_meas)+'\n'+'initial CM measurement mean:  '+str(CM_mu)+"\n"+'CM sigma (proportion):        '+str(CM_sigma_frac)+'\n'+'CM sigma (value):             '+str(CM_sigma_val)+'\n'+'\n'+'initial number of entities:   '+str(in_entities)+'\n'+'maximum number of entities:   '+str(max_entities)+'\n'+'\n'+'probability of fission:       '+str(prob_fis)+"\n"+'probability of sampling:      '+str(prob_dep)+"\n"+'\n')
	##q.close()
	
	# initialise array for saving output to file (+3 for n, x and y values)
	final_output = zeros((1, num_CM_meas+3))
	
	
	##################
	### INITIALISE ###
	##################
	
	# initial (empty) array for storing population density values in
	entities = in_entities
	
	# initial positions
	x = numpy.random.uniform(x_min, x_max, (entities,1))	
	y = numpy.random.uniform(y_min, y_max, (entities,1))	
	
	x_original = x.copy()
	y_original = y.copy()
	
	#plot(x,y,'.k')
	#axis([x_min, x_max, y_min, y_max])
	#hold(False)
	#show()
	
	# initial CM measurements
	CM_matrix = CM_mu*ones((entities,num_CM_meas))
	CM_sigma = CM_sigma_val*ones((num_CM_meas, ))
	
	
	###########################
	### START OF SIMULATION ###
	###########################
	
	start_time_2 = time.time()
	
	for n in range(generations):
		print 'generation: ', n, '\t', 'number of entities: ', entities
		
		### walk iteration ###
		xmove = numpy.random.normal(mig_mu, mig_sigma, (entities,1))	# amount to move by in x direction
		ymove = numpy.random.normal(mig_mu, mig_sigma, (entities,1))	# amount to move by in y direction
		
		# calculate new positions and check viability (geographical)
		xnew = x + xmove												# proposed x position
		xpos = where(less(xnew,x_min)|greater(xnew,x_max), x, xnew)		# check geographical viability of proposed x position
		x = xpos														# new x position
		
		ynew = y + ymove												# proposed  y position
		ypos = where(less(ynew,y_min)|greater(ynew,y_max), y, ynew)		# check geographical viability of proposed y position
		y = ypos														# new y position
		
		
		### fission / extinction iteration ###
		# entities undergo fission with probability == prob_fis
		CM_matrix_original = CM_matrix.copy()
		fis = scipy.random.random((entities,))												# entities undergo fission if value in corresponding fis array <= prob_fis
		
		# FISSION PROCESS - duplicate data for entities that undergo fission in x, y and CM_matrix
		fis_entities_ind = where(less_equal(fis,prob_fis))									# indices of entities undergoing fission
	#	print 'fis entities:', len(fis_entities_ind[0])
		
		x = vstack((x, x[fis_entities_ind]))
		y = vstack((y, y[fis_entities_ind]))
		CM_matrix = vstack((CM_matrix, CM_matrix_original[fis_entities_ind]))
		entities = x.shape[0]
		
		# check if number of entities has reached max_entities; if so:
		# EXTINCTION PROCESS - delete data for entities that undergo extinction from x, y and CM_matrix
	#	if entities>max_entities: print 'THERE ARE MORE ENTITIES THAN MAX ENTITIES - AN EXTINCTION PROCESS SHOULD BE HAPPENING HERE!'
	
		CM_matrix_original = CM_matrix.copy()
		
		if entities>max_entities:	
			ext_entities_ind = random.sample(range(entities), entities-max_entities)		# indices of surplus # of entities selected (random sampling) to undergo extinction
	#		if len(ext_entities_ind): print 'ext entities:', len(ext_entities_ind)
			
			x = delete(x, ext_entities_ind, axis=0)
			y = delete(y, ext_entities_ind, axis=0)
			CM_matrix = delete(CM_matrix_original, ext_entities_ind, axis=0)
			entities = x.shape[0]
		
		# check if all entities extinct
		if entities == 0:
			q = open('sim_data_2D_mig{0}_cm{1}_fis{2}_{3}.txt'.format("%.6f"% round(mig_sigma,12),CM_sigma_frac,prob_fis),taskID,'a')
			q.write("All entities have gone extinct at iteration "+str(n)+'.'+'\n')
			q.close()
			break
		
		
		### mutation iteration ###
		# NOTE! each CM measurement of each entity mutates (varies slightly) between generations
		# CM measurements vary according to gaussian distribution with mean CM_matrix and s.d. CM_sigma
		CM_matrix = scipy.random.normal(CM_matrix, CM_sigma, (entities,num_CM_meas))
		CM_matrix = abs(CM_matrix)
		
		
		### feature depositing iteration ###
		# entities sampled with probability == prob_dep
		# for sampled entities, data recoded to output: 25*n, x, y, CM_matrix values
		if (n>10000):
			dep = scipy.random.binomial(1, prob_dep, (entities,))
			dep_entities_ind = where(dep==1)												# indices of entities sampled
			dep_matrix = hstack((25*n*ones((len(dep_entities_ind[0]),1)), x[dep_entities_ind], y[dep_entities_ind], CM_matrix[dep_entities_ind]))		# 25*n, x, y and CM_matrix values for entities sampled
			final_output = vstack((final_output, dep_matrix))
	
	end_time_2 = time.time()
	time_taken_2 = end_time_2 - start_time_2
	print "time taken for", n, "generations is", time_taken_2/60, "minutes"
	
	#########################
	### END OF SIMULATION ###
	#########################
	
	
	##############
	### OUTPUT ###
	##############
	
	# record output to file
	q = open('sim_data_2D_mig{0}_cm{1}_fis{2}_{3}.txt'.format("%.6f"% round(mig_sigma,12),CM_sigma_frac,prob_fis,taskID),'a')
	savetxt(q, final_output[1:])
	q.close()
	
	run_r_simdata = 'R CMD BATCH "--args ' + str("%.6f"% round(mig_sigma,12) + '-' + str(CM_sigma_frac) + '-' + str(prob_fis) + '-' + str(taskID)) + '-' + str(generations) + '" simdata_2Dmodel.R'
	os.system(run_r_simdata)
	