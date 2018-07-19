# To run this, first, run source .bashrc & source .bash_profile to load up neuron library in python

from neuron import hoc,h
from netpyne import specs,sim

netParams_d1 = specs.NetParams()   # object of class NetParams to store the network parameters
netParams_d2 = specs.NetParams()   # object of class NetParams to store the network parameters

simConfig = specs.SimConfig()   # object of class SimConfig to store the simulation configuration


###############################################################################
# NETWORK PARAMETERS
#cellRule = netParams.importCellParams(label='PT5B_full', conds={'cellType': 'PT', 'cellModel': 'HH_full'},
#	  fileName='cells/PTcell.hoc', cellName='PTcell', cellArgs=[ihMod2str[cfg.ihModel], cfg.ihSlope], soma_0AtOrigin=True)
#	nonSpiny = ['apic_0', 'apic_1']
#netParams.popParams['PT5B'] = {'cellModel': 'HH_full', 'cellType': 'PT', 'ynormRange': layer['5B'], 'numCells':1}
###############################################################################

# Population parameters
cellRule_d1=netParams_d1.importCellParams(label='D1MSN', conds={'cellType': 'D1', 'cellModel': 'D1MSN'}, fileName='d1msn.hoc', cellName='d1msn')
cellRule_d2=netParams_d2.importCellParams(label='D2MSN', conds={'cellType': 'D2', 'cellModel': 'D2MSN'}, fileName='d2msn.hoc', cellName='d2msn')

netParams_d1.popParams['D1MSN'] = {'cellModel': 'D1MSN', 'cellType': 'D1', 'numCells': 10} # add dict with params for this pop 
netParams_d2.popParams['D2MSN'] = {'cellModel': 'D2MSN', 'cellType': 'D2', 'numCells': 1} # add dict with params for this pop 

netParams_d1.stimSourceParams['AMPA'] = {'type':'NetStim','rate':10,'noise':0.5}
netParams_d1.stimTargetParams['AMPA->D1MSN'] = {'source': 'AMPA', 'conds': {'cellType': 'D1'}, 'weight': 0.01, 'delay': 5, 'synMech': 'exc'}

##################
# Cell parameters
#cellRule_d1 = {'conds': {'cellModel': 'D1MSN', 'cellType': 'D1'},  'secs': {}} 	# cell rule dict
#cellRule_d2 = {'conds': {'cellModel': 'd2msn', 'cellType': 'D2MSN'},  'secs': {}} 	# cell rule dict
# This is how to manipulate cellular parameters but we do not need these now
#cellRule_d1['secs']['soma_0'] = {'geom': {}, 'mechs': {}}  														# soma_0 params dict
#cellRule_d1['secs']['soma_0']['geom'] = {'diam': 18.8, 'L': 18.8, 'Ra': 123.0}  									# soma_0 geometry
#cellRule['secs']['soma_0']['mechs']['d2msn'] = {'gnabar': 0.12, 'gkbar': 0.036, 'gl': 0.003, 'el': -70}  		# soma_0 hh mechanism
###################
cellRule_d1['secs']['soma_0']['vinit'] = -80
cellRule_d1['secs']['soma_0']['vinit'] = -80
netParams_d1.cellParams['D1MSN'] = cellRule_d1  												# add dict to list of cell params
netParams_d2.cellParams['D2MSN'] = cellRule_d2 

# Synaptic mechanism parameters
netParams_d1.synMechParams['AMPA'] = {'mod': 'AMPA', 'tau_r': 15.3604, 'tau_d': 1.55121, 'gbar': 8.5e-04}
netParams_d1.synMechParams['GABA'] = {'mod': 'GABA', 'tau_r': 2.3578, 'tau_d': 7.7747}

netParams_d2.synMechParams['AMPA'] = {'mod': 'AMPA', 'tau_r': 15.3604, 'tau_d': 1.55121, 'gbar': 9e-04}
netParams_d2.synMechParams['GABA'] = {'mod': 'GABA', 'tau_r': 2.3578, 'tau_d': 7.7747}

s
# Stimulation parameters
netParams_d1.stimSourceParams['AMPA'] = {'type': 'NetStim', 'rate': 10, 'noise': 0.5, 'start': 1}
netParams_d1.stimTargetParams['AMPA->D1MSN'] = {'source': 'AMPA', 'conds': {'pop': 'D1MSN'}, 'weight': 0.1, 'delay': 'uniform(1,5)'}

netParams_d2.stimSourceParams['AMPA'] = {'type': 'NetStim', 'rate': 10, 'noise': 0.5, 'start': 1}
netParams_d2.stimTargetParams['AMPA->D2MSN'] = {'source': 'AMPA', 'conds': {'pop': 'D2MSN'}, 'weight': 0.1, 'delay': 'uniform(1,5)'}

cellRule_d1['secs']['soma_0']['spikeGenLoc'] = 1.0
cellRule_d2['secs']['soma_0']['spikeGenLoc'] = 1.0

netParams_d1.stimSourceParams['Input_1'] = {'type': 'IClamp', 'del': 300, 'dur': 100, 'amp': 200}
netParams_d1.stimSourceParams['Input_4'] = {'type': 'NetStim', 'interval': 'uniform(20,100)', 'start': 600, 'noise': 0.1}
netParams_d1.stimTargetParams['Input_1->D1MSN'] = {'source': 'Input_1', 'sec':'soma_0', 'loc': 0.8, 'conds': {'pop':'D1MSN', 'cellList': range(10)}}
netParams_d1.stimTargetParams['Input_4->D1MSN'] = {'source': 'Input_4', 'sec':'soma_0', 'loc': 0.5, 'weight': '0.1+normal(0.2,0.05)','delay': 1,
                              				'conds': {'cellType':'D1', 'ynorm': [0.6,1.0]}}
netParams_d1.stimSourceParams['bkg'] = {'type': 'NetStim', 'rate': 10, 'noise': 0.5}
netParams_d1.stimTargetParams['bkg->D1'] = {'source': 'bkg', 'conds': {'cellType': 'D1'}, 'weight': 0.01, 'delay': 5, 'synMech': 'exc'}

# Connectivity parameters
netParams_d1.connParams['D1MSN->D1MSN'] = {
    'preConds': {'pop': 'D1MSN'}, 
    'postConds': {'pop': 'D1MSN'},
    'sec':'soma_0',
    'synMech':'AMPA',
    'weight': 0.01,                    # weight of each connection
    'delay': '0.1+normal(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
    'probability': 0.4}
    
netParams_d2.connParams['D1MSN->D2MSN'] = {
    'preConds': {'pop': 'D1MSN'}, 
    'postConds': {'pop': 'D2MSN'},
    'sec':'soma_0',
    'synMech':'GABA',
    'weight': 0.01,                    # weight of each connection
    'delay': '0.2+normal(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
    'probability': 0.4}

netParams_d2.connParams['D2MSN->D1MSN'] = {
    'preConds': {'pop': 'D2MSN'}, 
    'postConds': {'pop': 'D1MSN'},
    'sec':'proximal',
    'synMech':'AMPA',
    'weight': 0.01,                    # weight of each connection
    'delay': '0.2+normal(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
    'probability': 0.4}
    
netParams_d2.connParams['D2MSN->D1MSN'] = {
    'preConds': {'pop': 'D2MSN'}, 
    'postConds': {'pop': 'D1MSN'},
    'sec':'proximal',
    'synMech':'GABA',
    'weight': 0.01,                    # weight of each connection
    'delay': '0.2+normal(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
    'probability': 0.4}

##########
#netParams_d2.connParams['D2MSN->D2MSN'] = {
#    'preConds': {'pop': 'D2MSN'}, 
#    'postConds': {'pop': 'D2MSN'},
#    'sec':'dend',
#    'synMech':'AMPA',
#    'weight': 0.01,                    # weight of each connection
#    'delay': '0.2+normal(13.0,1.4)',     # delay min=0.2, mean=13.0, var = 1.4
#    'probability': 0.4}                    # probability of connection
##########

###############################################################################
# SIMULATION PARAMETERS
###############################################################################

# Simulation parameters
simConfig.duration = 1*200 # Duration of the simulation, in ms
simConfig.dt = 0.025 # Internal integration timestep to use
simConfig.seeds = {'conn': 1, 'stim': 1, 'loc': 1} # Seeds for randomizers (connectivity, input stimulation and cell locations)
simConfig.createNEURONObj = 1  # create HOC objects when instantiating network
simConfig.createPyStruct = 1  # create Python structure (simulator-independent) when instantiating network
simConfig.verbose = True  # show detailed messages 


# Recording 
simConfig.recordCells = [10]  # which cells to record from
simConfig.recordTraces = {'Vsoma_0':{'sec':'soma_0','loc':0.5,'var':'v'}}
simConfig.recordStim = True  # record spikes of cell stims
simConfig.recordStep = 0.1 # Step size in ms to save data (eg. V traces, LFP, etc)

# Saving
simConfig.filename = 'msn_net'  # Set file output name
simConfig.saveFileStep = 1000 # step size in ms to save data to disk
simConfig.savePickle = False # Whether or not to write spikes etc. to a .mat file

# Analysis and plotting 
simConfig.analysis['plotRaster'] = True  # Plot raster
simConfig.analysis['plotTraces'] = {'include': [('D1MSN',0)]}  
# simConfig.analysis['plotTraces'] = {'include': [('D1MSN',0)]}
simConfig.analysis['plot2Dnet'] = True  # Plot 2D net cells and connections

#simConfig.recordLFP = [[-15, y, 1.0*netParams_d2.sizeZ] for y in range(netParams_d2.sizeY/5, netParams_d2.sizeY, netParams_d2.sizeY/5)]
#simConfig.analysis['plotLFP'] = True

sim.createSimulateAnalyze(netParams_d1, simConfig)
