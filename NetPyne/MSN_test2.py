
from neuron import hoc,h
from netpyne import specs,sim
from netpyne import utils

netParams_d2 = specs.NetParams()   # object of class NetParams to store the network parameters
netParams_d2.popParams['D2MSN'] = {'cellModel': 'MSN', 'cellType': 'D2', 'numCells': 1} # add dict with params for this pop 

cellRule_d2=netParams_d2.importCellParams(label='D2MSN', conds={'cellType': 'D2', 'cellModel': 'MSN'}, fileName='d2msn.py', cellName='WTD2')


netParams_d2.cellParams['D2MSN'] = cellRule_d2  

#utils.mechVarList()
# check different mechs
# netParams_d2.cellParams['D2MSN']['secs']['soma_0']['mechs']['naf']
#for secName in cellRule_d2['secs']:
#	cellRule_d2['secs'][secName]['mechs']['kaf'].qfact=1.2

netParams_d2.stimSourceParams['Input_1'] = {'type': 'IClamp', 'del': 100, 'dur': 2000, 'amp': 0.08}
netParams_d2.stimTargetParams['Input_1->D2MSN'] = {'source': 'Input_1','sec':'soma_0', 'loc': 0.5, 'conds': {'pop':'D2MSN', 'cellList': [0]}}

simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration

simConfig.duration = 2515        # Duration of the simulation, in ms
simConfig.dt = 0.25                # Internal integration timestep to use
simConfig.verbose = True           # Show detailed messages 
simConfig.recordTraces = {'v_soma':{'sec':'soma_0','loc': 0.5, 'var':'v'}}
simConfig.recordStep = 1         # Step size in ms to save data (eg. V traces, LFP, etc)
simConfig.filename = 'model_output'  # Set file output name
simConfig.savePickle = False      # Save params, network and sim output to pickle file
simConfig.saveJson = True
#,'hshift_caL':0,'hqfact_caL13':1.2,'hshift_caL13':0
simConfig.hParams = {'v_init': -80.0}
simConfig.analysis['plotTraces'] = {'include': [0]} 			# Plot recorded traces for this list of cells
# simConfig.analysis['plotShape'] = True
sim.createSimulateAnalyze(netParams = netParams_d2, simConfig = simConfig)    

#for sec in h.allsec():
#	print sec.name(),sec.nseg
#
#for sec in h.soma():
#	print sec.name(),sec.ek