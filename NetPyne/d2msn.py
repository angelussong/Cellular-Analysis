import sys
from neuron import h, hoc
from math import sqrt, pi, log, exp
Ra=100
ek=-94
ena=50
with open('wtd2_params.txt') as f:
	for line in f:
		exec(line)

h.xopen("msn_wolf_all_tau_vecs.hoc")
dgnaf_proxF = G_NAFD / G_NAF
dgnap_proxF = G_NAPD / G_NAP
dgnaf_midF = G_NAFD / G_NAF
dgnap_midF = G_NAPD / G_NAP
dgkas_midF = G_KASD / G_KAS
dgkaf_midF = G_KAFD / G_KAF
dgnaf_distF = G_NAFD / G_NAF
dgnap_distF = G_NAPD / G_NAP
dgkas_distF = G_KASD / G_KAS
dgkaf_distF = G_KAFD / G_KAF

class WTD2:
	def __init__(self):
		self.create_cell()

		self.optimize_nseg()
		self.settopo()
		self.add_all()
		self.geom_nseg()
		self.setchannels()

	
	def geom_nseg (self):
		freq = 100 #Hz
		d_lambda = 0.1
		for sec in self.all: 
			sec.nseg = int((sec.L/(d_lambda*self.lambda_f(sec))+0.9)/2)*2 + 1
		for sec in self.dist:
			sec.nseg=3
		for sec in self.axon:
			sec.nseg=1
		for sec in self.myelin:
			sec.nseg=1 
		
	def lambda_f (self, section):
		# these are reasonable values for most models
		freq = 100
		# The lowest number of n3d() is 2
		if (section.n3d() < 2):
			return 1e5*sqrt(section.diam/(4*pi*freq*section.Ra*section.cm))
		# above was too inaccurate with large variation in 3d diameter
		# so now we use all 3-d points to get a better approximate lambda
		x1 = section.arc3d(0)
		d1 = section.diam3d(0)
		self.lam = 0
		#print section, " n3d:", section.n3d(), " diam3d:", section.diam3d(0)
		for i in range(section.n3d()): #h.n3d()-1
			x2 = section.arc3d(i)
			d2 = section.diam3d(i)
			self.lam += (x2 - x1)/sqrt(d1 + d2)
			x1 = x2
			d1 = d2
			#  length of the section in units of lambda
		self.lam *= sqrt(2) * 1e-5*sqrt(4*pi*freq*section.Ra*section.cm)
		return section.L/self.lam

	def create_cell(self):
		self.soma = [h.Section(name='soma[%d]' % i) for i in xrange(1)]
		self.proximal = [h.Section(name='proximal[%d]' % i) for i in xrange(6)]
		self.middend = [h.Section(name='middend[%d]' % i) for i in xrange(12)]
		self.distal = [h.Section(name='distal[%d]' % i) for i in xrange(26)]
		self.axon = [h.Section(name='axon[%d]' % i) for i in xrange(2)]
		self.myelin = [h.Section(name='myelin[%d]' % i) for i in xrange(1)] 
		
		h.pt3dclear(sec = self.soma[0])
		h.pt3dadd(0,0,0,1, sec = self.soma[0])
		h.pt3dadd(15,0,0,1, sec = self.soma[0])
		
		self.proximal[0].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[0])
		h.pt3dadd(15,0,0,1, sec=self.proximal[0])
		h.pt3dadd(45,-29,0,1, sec=self.proximal[0])
		
		self.proximal[1].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[1])
		h.pt3dadd(15,0,0,1, sec=self.proximal[1])
		h.pt3dadd(45,30,0,1, sec=self.proximal[1])
		
		self.proximal[2].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[2])
		h.pt3dadd(0,0,0,1, sec=self.proximal[2])
		h.pt3dadd(-29,-29,0,1, sec=self.proximal[2])
		
		self.proximal[3].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[3])
		h.pt3dadd(0,0,0,1, sec=self.proximal[3])
		h.pt3dadd(-29,30,0,1, sec=self.proximal[3])
		
		self.proximal[4].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[4])
		h.pt3dadd(15,0,0,1, sec=self.proximal[4])
		h.pt3dadd(45,0,0,1, sec=self.proximal[4])
		
		self.proximal[5].connect(self.soma[0])
		h.pt3dclear(sec = self.proximal[5])
		h.pt3dadd(0,0,0,1, sec=self.proximal[5])
		h.pt3dadd(-30,0,0,1, sec=self.proximal[5])
		
		# Finished with proximal Now to middend
		
		self.middend[0].connect(self.proximal[0])
		h.pt3dclear(sec = self.middend[0])
		h.pt3dadd(45,-29,0,1, sec=self.middend[0])
		h.pt3dadd(45,-59,0,1, sec=self.middend[0])
		
		self.middend[1].connect(self.proximal[0])
		h.pt3dclear(sec = self.middend[1])
		h.pt3dadd(45,-29,0,1, sec=self.middend[1])
		h.pt3dadd(90,-29,0,1, sec=self.middend[1])
		
		self.middend[2].connect(self.proximal[1])
		h.pt3dclear(sec = self.middend[2])
		h.pt3dadd(45,30,0,1, sec=self.middend[2])
		h.pt3dadd(90,30,0,1, sec=self.middend[2])
		
		self.middend[3].connect(self.proximal[1])
		h.pt3dclear(sec = self.middend[3])
		h.pt3dadd(45,30,0,1, sec=self.middend[3])
		h.pt3dadd(60,75,0,1, sec=self.middend[3])
		
		self.middend[4].connect(self.proximal[2])
		h.pt3dclear(sec = self.middend[4])
		h.pt3dadd(-29,-29,0,1, sec=self.middend[4])
		h.pt3dadd(-44,-59,0,1, sec=self.middend[4])
		
		self.middend[5].connect(self.proximal[2])
		h.pt3dclear(sec = self.middend[5])
		h.pt3dadd(-29,-29,0,1, sec=self.middend[5])
		h.pt3dadd(-74,-29,0,1, sec=self.middend[5])
		
		self.middend[6].connect(self.proximal[3])
		h.pt3dclear(sec = self.middend[6])
		h.pt3dadd(-29,30,0,1, sec=self.middend[6])
		h.pt3dadd(-44,75,0,1, sec=self.middend[6])
		
		self.middend[7].connect(self.proximal[3])
		h.pt3dclear(sec = self.middend[7])
		h.pt3dadd(-29,30,0,1, sec=self.middend[7])
		h.pt3dadd(-74,30,0,1, sec=self.middend[7])
		
		self.middend[8].connect(self.proximal[4])
		h.pt3dclear(sec = self.middend[8])
		h.pt3dadd(45,  0, 0, 1, sec=self.middend[8])
		h.pt3dadd(75, 15, 0, 1, sec=self.middend[8])
		
		self.middend[9].connect(self.proximal[4])
		h.pt3dclear(sec = self.middend[9])
		h.pt3dadd(45,  0, 0, 1, sec=self.middend[9])
		h.pt3dadd(75, -15, 0, 1, sec=self.middend[9])
		
		self.middend[10].connect(self.proximal[5])
		h.pt3dclear(sec = self.middend[10])
		h.pt3dadd(-30, 0, 0, 1, sec=self.middend[10])
		h.pt3dadd(-60, 15, 0, 1, sec=self.middend[10])
		
		self.middend[11].connect(self.proximal[5])
		h.pt3dclear(sec = self.middend[11])
		h.pt3dadd(-30, 0, 0, 1, sec=self.middend[11])
		h.pt3dadd(-60, -15, 0, 1, sec=self.middend[11])
		
		# Finished with middend Now to distal
		self.distal[0].connect(self.middend[0])
		h.pt3dclear(sec=self.distal[0])
		h.pt3dadd(45, -59, 0, 1, sec=self.distal[0])
		h.pt3dadd(30, -89, 0, 1, sec=self.distal[0])
		
		self.distal[1].connect(self.middend[0])
		h.pt3dclear(sec=self.distal[1])
		h.pt3dadd(45, -59, 0, 1, sec=self.distal[1])
		h.pt3dadd(75, -74, 0, 1, sec=self.distal[1])
		
		self.distal[2].connect(self.middend[1])
		h.pt3dclear(sec=self.distal[2])
		h.pt3dadd(90, -29, 0, 1, sec=self.distal[2])
		h.pt3dadd(120, -14, 0, 1, sec=self.distal[2])
		
		self.distal[3].connect(self.middend[1])
		h.pt3dclear(sec=self.distal[3])
		h.pt3dadd(90, -29, 0, 1, sec=self.distal[3])
		h.pt3dadd(120, -59, 0, 1, sec=self.distal[3])
		
		self.distal[4].connect(self.middend[2])
		h.pt3dclear(sec=self.distal[4])
		h.pt3dadd(90, 30, 0, 1, sec=self.distal[4])
		h.pt3dadd(120, 15, 0, 1, sec=self.distal[4])
		
		self.distal[5].connect(self.middend[2])
		h.pt3dclear(sec=self.distal[5])
		h.pt3dadd(90, 30, 0, 1, sec=self.distal[5])
		h.pt3dadd(105, 60, 0, 1, sec=self.distal[5])
		
		self.distal[6].connect(self.middend[3])
		h.pt3dclear(sec=self.distal[6])
		h.pt3dadd(60, 75, 0, 1, sec=self.distal[6])
		h.pt3dadd(90, 90, 0, 1, sec=self.distal[6])
		
		self.distal[7].connect(self.middend[3])
		h.pt3dclear(sec=self.distal[7])
		h.pt3dadd(60, 75, 0, 1, sec=self.distal[7])
		h.pt3dadd(45, 105, 0, 1, sec=self.distal[7])
		
		self.distal[8].connect(self.middend[4])
		h.pt3dclear(sec=self.distal[8])
		h.pt3dadd(-44, -59, 0, 1, sec=self.distal[8])
		h.pt3dadd(-29, -89, 0, 1, sec=self.distal[8])
		
		self.distal[9].connect(self.middend[4])
		h.pt3dclear(sec=self.distal[9])
		h.pt3dadd(-44, -59, 0, 1, sec=self.distal[9])
		h.pt3dadd(-89, -74, 0, 1, sec=self.distal[9])
		
		self.distal[10].connect(self.middend[5])
		h.pt3dclear(sec=self.distal[10])
		h.pt3dadd(-74, -29, 0, 1, sec=self.distal[10])
		h.pt3dadd(-104, -14, 0, 1, sec=self.distal[10])
		
		self.distal[11].connect(self.middend[5])
		h.pt3dclear(sec=self.distal[11])
		h.pt3dadd(-74, -29, 0, 1, sec=self.distal[11])
		h.pt3dadd(-104, -59, 0, 1, sec=self.distal[11])
		
		self.distal[12].connect(self.middend[6])
		h.pt3dclear(sec=self.distal[12])
		h.pt3dadd(-44, 75, 0, 1, sec=self.distal[12])
		h.pt3dadd(-14, 105, 0, 1, sec=self.distal[12])
		
		self.distal[13].connect(self.middend[6])
		h.pt3dclear(sec=self.distal[13])
		h.pt3dadd(-44, 75, 0, 1, sec=self.distal[13])
		h.pt3dadd(-74, 90, 0, 1, sec=self.distal[13])
		
		self.distal[14].connect(self.middend[7])
		h.pt3dclear(sec=self.distal[14])
		h.pt3dadd(-74, 30, 0, 1, sec=self.distal[14])
		h.pt3dadd(-89, 60, 0, 1, sec=self.distal[14])
		
		self.distal[15].connect(self.middend[7])
		h.pt3dclear(sec=self.distal[15])
		h.pt3dadd(-74, 30, 0, 1, sec=self.distal[15])
		h.pt3dadd(-104, 15, 0, 1, sec=self.distal[15])
		
		self.distal[16].connect(self.middend[8])
		h.pt3dclear(sec=self.distal[16])
		h.pt3dadd(75, 15, 0, 1, sec=self.distal[16])
		h.pt3dadd(90, 15, 0, 1, sec=self.distal[16])
		
		self.distal[17].connect(self.middend[8])
		h.pt3dclear(sec=self.distal[17])
		h.pt3dadd(75, 15, 0, 1, sec=self.distal[17])
		h.pt3dadd(90, 30, 0, 1, sec=self.distal[17])
		
		self.distal[18].connect(self.middend[9])
		h.pt3dclear(sec=self.distal[18])
		h.pt3dadd(75, -15, 0, 1, sec=self.distal[18])
		h.pt3dadd(90, -15, 0, 1, sec=self.distal[18])
		
		self.distal[19].connect(self.middend[9])
		h.pt3dclear(sec=self.distal[19])
		h.pt3dadd(75, -15, 0, 1, sec=self.distal[19])
		h.pt3dadd(90, -30, 0, 1, sec=self.distal[19])
		
		self.distal[20].connect(self.middend[10])
		h.pt3dclear(sec=self.distal[20])
		h.pt3dadd(-40, 15, 0, 1, sec=self.distal[20])
		h.pt3dadd(-63, 18, 0, 1, sec=self.distal[20])
		
		self.distal[21].connect(self.middend[10])
		h.pt3dclear(sec=self.distal[21])
		h.pt3dadd(-60, 15, 0, 1, sec=self.distal[21])
		h.pt3dadd(-75, 30, 0, 1, sec=self.distal[21])
		# this is problematic! 
		self.distal[22].connect(self.middend[10],0.8,0)
		h.pt3dclear(sec=self.distal[22])
		h.pt3dadd(-54, 12, 0, 1, sec=self.distal[22])
		h.pt3dadd(-69, 27, 0, 1, sec=self.distal[22])
		
		self.distal[23].connect(self.middend[11])
		h.pt3dclear(sec=self.distal[23])
		h.pt3dadd(-40, 0, 0, 1, sec=self.distal[23])
		h.pt3dadd(-63, -12, 0, 1, sec=self.distal[23])
		
		self.distal[24].connect(self.middend[11])
		h.pt3dclear(sec=self.distal[24])
		h.pt3dadd(-40, 0, 0, 1, sec=self.distal[24])
		h.pt3dadd(-63, -18, 0, 1, sec=self.distal[24])
		# this is problematic
		self.distal[25].connect(self.middend[11],0.8,0)
		h.pt3dclear(sec=self.distal[25])
		h.pt3dadd(-54, -12, 0, 1, sec=self.distal[25])
		h.pt3dadd(-69, -27, 0, 1, sec=self.distal[25])
		
		self.axon[0].connect(self.soma[0])
		h.pt3dclear(sec=self.axon[0])
		h.pt3dadd(7.5,0,0,1, sec=self.axon[0])
		h.pt3dadd(7.5,20,0,1, sec=self.axon[0])
		
		self.axon[1].connect(self.axon[0])
		h.pt3dclear(sec=self.axon[1])
		h.pt3dadd(7.5,20,0,1, sec=self.axon[1])
		h.pt3dadd(7.5,50,0,1, sec=self.axon[1])
		
		self.myelin[0].connect(self.axon[1])
		h.pt3dclear(sec=self.myelin[0])
		h.pt3dadd(7.5,50,0,1, sec=self.myelin[0])
		h.pt3dadd(7.5,450,0,1, sec=self.myelin[0])


	def optimize_nseg (self):
		self.all=h.SectionList()
		self.prox=h.SectionList()
		self.mid=h.SectionList()
		self.dist=h.SectionList()
		self.all_no_axon=h.SectionList()
		for section in self.soma:
			self.all.append(sec=section)
			self.all_no_axon.append(sec=section)
		for section in self.proximal:
			self.all.append(sec=section)
			self.all_no_axon.append(sec=section)
			self.prox.append(sec=section)
		for section in self.middend:
			self.all.append(sec=section)
			self.all_no_axon.append(sec=section)
			self.mid.append(sec=section)
		for section in self.distal:
			self.all.append(sec=section)
			self.all_no_axon.append(sec=section)
			self.dist.append(sec=section)
		for section in self.axon:
			self.all.append(sec=section)
		for section in self.myelin:
			self.all.append(sec=section)
		
		for sec in self.all:
			sec.Ra=100
			sec.cm=0.900001
		for sec in self.myelin:
			sec.cm=0.01800002
		
	
	def settopo(self):
		#### set L and diam
		self.soma[0].L=10.4941
		self.soma[0].diam=9.5896

		for sec in self.proximal:
			sec.L=23.7905
			sec.diam=1.2933559
		for sec in self.middend:
			sec.L=26.846188
			sec.diam=1.2744389
		for sec in self.distal:
			sec.L=116.49579
			sec.diam=0.90830576
		self.axon[0].L=20
		self.axon[0].diam=1.1875
		self.axon[1].L=30
		self.axon[1].diam=0.6875
		self.myelin[0].L=400
		self.myelin[0].diam=0.5
  
	def add_all(self):
		for sec in self.all:
			sec.insert('pas')
		for sec in self.all_no_axon:
			sec.insert('naf')
			sec.insert('nap')
			sec.insert('kir')
			sec.insert('kas')
			sec.insert('kaf')
			sec.insert('bkkca')
			sec.insert('skkca')
			sec.insert('kdr')
			sec.insert('caldyn')
			sec.insert('caL')
			sec.insert('caL13')
			sec.insert('cadyn')
			sec.insert('can')
			sec.insert('car')
			sec.insert('cat')
		for sec in self.axon:
			sec.insert('naf')
			sec.insert('nap')
			sec.insert('kas')
			sec.insert('kaf')
			sec.insert('kdr')
		
	def setchannels(self):
		#### set pas params
		for sec in self.all:
			sec.g_pas=5.65446e-05
			sec.e_pas=-76.944364
		#### set kir params
		for sec in self.all_no_axon:
			sec.qfact_kir=0.2
			sec.mshift_kir=SHKIR
			sec.mslope_kir=11.889457
			sec.gkbar_kir=0.013600679
			sec.cai=5.5097902e-06
			sec.cao=5
			sec.cali=5.5139085e-06
			sec.calo=5
			sec.eca=165.15653
			sec.ena=50
			sec.ek=-94
		#### set shift and dynamic params
		for sec in self.all_no_axon:
			sec.mshift_kas = MSHIFT_KAS
			sec.hshift_kas = HSHIFT_KAS
			sec.mshift_kaf = MSHIFT_KAF
			sec.hshift_kaf = HSHIFT_KAF
		for sec in self.axon:
			sec.mshift_naf=-10
			sec.hshift_naf=-5
			sec.ena=50
			sec.ek=-94
		for sec in self.all_no_axon:
			QSCALE=0.2/0.5
			sec.mqfact_naf=3*QSCALE
			sec.hqfact_naf=3*QSCALE
			sec.qfact_kas=9*QSCALE
			sec.qfact_kaf=3*QSCALE
			sec.qfact_nap=0.5
			sec.qfact_kdr=0.5*QSCALE
			sec.qfact_skkca = 1 * QSCALE
			sec.q_bkkca = 1 * QSCALE
			sec.qfact_can = 3 * QSCALE
			sec.qfact_car = 3 * QSCALE
			sec.qfact_cat = 3 * QSCALE
			sec.qfact_caL = 3 * QSCALE
			sec.qfact_caL13 = 3 * QSCALE
			sec.hqfact_caL13 = 3 * QSCALE         
		for sec in self.axon:
			QSCALE=0.2/0.5
			sec.mqfact_naf=3*QSCALE
			sec.hqfact_naf=3*QSCALE
			sec.qfact_kas=9*QSCALE
			sec.qfact_kaf=3*QSCALE
			sec.qfact_nap=0.5
			sec.qfact_kdr=0.5*QSCALE
		##### set conductance in soma
		for sec in self.soma:
			# sec.gnabar_naf=0.890840924
			sec.gnabar_naf=0.900840924
			# sec.gnabar_naf=0
			# sec.gnabar_nap=0
			sec.gnabar_nap=0.0004171898
			sec.gkbar_kas=0.0000954872
			sec.gkbar_kaf=0.0306049324
			# sec.gkbar_kas=0
			# sec.gkbar_kaf=0
			sec.gkbar_kdr=0.0000395576
			# sec.gkbar_kdr=0
			sec.gkbar_bkkca=0.0005419430
			sec.gkbar_skkca=0.0634955575
			sec.pbar_caL=0.0000019430
			sec.pcaLbar_caL13=0.0000001262
			sec.pbar_can=0.0000037525
			sec.pcarbar_car=0.0000116228 
			sec.pcatbar_cat=0.0000001111
			# sec.gkbar_bkkca=0.
			# sec.gkbar_skkca=0
			# sec.pbar_caL=0
			# sec.pcaLbar_caL13=0
			# sec.pbar_can=0
			# sec.pcarbar_car=0
			# sec.pcatbar_cat=0
		for sec in self.prox:
			# sec.gnabar_naf=round(0.890840924*dgnaf_proxF,9)
			sec.gnabar_naf=0.900840924*dgnaf_proxF
			# sec.gnabar_naf=0
			# sec.gnabar_nap=0
			sec.gnabar_nap=0.0004171898*dgnap_proxF
			sec.gkbar_kas=0.0000954872
			sec.gkbar_kaf=0.0306049324
			# sec.gkbar_kas=0
			# sec.gkbar_kaf=0
			sec.gkbar_kdr=0.0000395576
			# sec.gkbar_kdr=0
			sec.gkbar_bkkca=0.0005419430
			sec.gkbar_skkca=0.0634955575
			sec.pbar_caL=0.0000019430
			sec.pcaLbar_caL13=0.0000001262
			sec.pbar_can=0.0000037525
			sec.pcarbar_car=0.0000116228 
			sec.pcatbar_cat=0.0000001111
			# sec.gkbar_bkkca=0.
			# sec.gkbar_skkca=0
			# sec.pbar_caL=0
			# sec.pcaLbar_caL13=0
			# sec.pbar_can=0
			# sec.pcarbar_car=0
			# sec.pcatbar_cat=0
		for sec in self.mid:
			# sec.gnabar_naf=round(0.890840924*dgnaf_midF,9)
			sec.gnabar_naf=0.900840924*dgnaf_midF
			# sec.gnabar_naf=0
			# sec.gnabar_nap=0
			sec.gnabar_nap=0.0004171898*dgnap_midF
			sec.gkbar_kas=0.0000954872*dgkas_midF
			sec.gkbar_kaf=0.0306049324*dgkaf_midF
			# sec.gkbar_kas=0
			# sec.gkbar_kaf=0
			sec.gkbar_kdr=0.0000395576
			# sec.gkbar_kdr=0
			sec.gkbar_bkkca=0.0005419430
			sec.gkbar_skkca=0.0634955575
			sec.pbar_caL=0.0000019430
			sec.pcaLbar_caL13=0.0000001262
			sec.pbar_can=0.0000037525
			sec.pcarbar_car=0.0000116228 
			sec.pcatbar_cat=0.0000001111
			# sec.gkbar_bkkca=0
			# sec.gkbar_skkca=0
			# sec.pbar_caL=0
			# sec.pcaLbar_caL13=0
			# sec.pbar_can=0
			# sec.pcarbar_car=0
			# sec.pcatbar_cat=0
		for sec in self.dist:
			# sec.gnabar_naf=0.890840924*dgnaf_distF
			sec.gnabar_naf=0.90840924*dgnaf_distF
			# sec.gnabar_naf=0
			# sec.gnabar_nap=0
			sec.gnabar_nap=0.0004171898*dgnap_distF
			sec.gkbar_kas=0.0000954872*dgkas_distF
			sec.gkbar_kaf=0.0306049324*dgkaf_distF
			sec.qfact_kaf=0.75
			# sec.gkbar_kas=0
			# sec.gkbar_kaf=0
			sec.gkbar_kdr=0.0000395576
			# sec.gkbar_kdr=0
			sec.gkbar_bkkca=0.0005419430
			sec.gkbar_skkca=0.0634955575
			sec.pbar_caL=0.0000019430
			sec.pcaLbar_caL13=0.0000001262
			sec.pbar_can=0.0000037525
			sec.pcarbar_car=0.0000116228 
			sec.pcatbar_cat=0.0000001111
			# sec.gkbar_bkkca=0.
			# sec.gkbar_skkca=0
			# sec.pbar_caL=0
			# sec.pcaLbar_caL13=0
			# sec.pbar_can=0
			# sec.pcarbar_car=0
			# sec.pcatbar_cat=0
		for sec in self.axon:
			sec.gnabar_naf=1.5652924106
			sec.gkbar_kas=0.0267535675
			sec.gkbar_kaf=0.0047530889
			sec.gnabar_nap=0.00004
			sec.gkbar_kdr=0.0001677798
			# sec.gnabar_naf=0
			# sec.gkbar_kas=0
			# sec.gkbar_kaf=0
			# sec.gnabar_nap=0
			# sec.gkbar_kdr=0



	    
	    
    	
    		
    	
    	
    	
        
        
        
        
        
        
        
