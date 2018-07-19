
from neuron import hoc,h
from netpyne import specs,sim
from netpyne import utils

netParams_d2 = specs.NetParams()   # object of class NetParams to store the network parameters

cellRule_d2=netParams_d2.importCellParams(label='D2MSN', conds={'cellType': 'D2', 'cellModel': 'MSN'}, fileName='d2msn_gai.hoc', cellName='d2msn')

netParams_d2.popParams['D2MSN'] = {'cellModel': 'MSN', 'cellType': 'D2', 'numCells': 1} # add dict with params for this pop 

netParams_d2.cellParams['D2MSN'] = cellRule_d2  
for secName in cellRule_d2['secs']:
 	cellRule_d2['secs'][secName]['mechs']['naf'] 

utils.mechVarList()
# check different mechs
# netParams_d2.cellParams['D2MSN']['secs']['soma']['mechs']['naf']

netParams_d2.stimSourceParams['Input_1'] = {'type': 'IClamp', 'del': 0, 'dur': 400, 'amp': 0}
netParams_d2.stimTargetParams['Input_1->D2MSN'] = {'source': 'Input_1','sec':'soma', 'loc': 0.5, 'conds': {'pop':'D2MSN', 'cellList': range(1)}}

simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration

simConfig.duration = 0.001       # Duration of the simulation, in ms
simConfig.dt = 0.00001                # Internal integration timestep to use
simConfig.verbose = True           # Show detailed messages 
simConfig.recordTraces = {'mnaf_soma':{'sec':'soma','loc': 0, 'var':'m_naf'},'hnaf_soma':{'sec':'soma','loc':0,'var':'h_naf'}}
simConfig.recordStep = 0.00001         # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.filename = 'model_output'  # Set file output name
simConfig.savePickle = False      # Save params, network and sim output to pickle file
simConfig.saveJson = True
#,'hshift_caL':0,'hqfact_caL13':1.2,'hshift_caL13':0
# simConfig.hParams = {'celsius': 35, 'v_init': -80.0}
simConfig.analysis['plotTraces'] = {'include': [0]} 			# Plot recorded traces for this list of cells

sim.createSimulateAnalyze(netParams = netParams_d2, simConfig = simConfig)    
