import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

#import microwave_networks as mw_net
import skmd as md


mpl.rc('axes',labelsize=24)
mpl.rc('font',size=24)

### Define constant
c = 3e8 # Velocity of light in meters/second

epsilon_r = 10 # dielectric constant or the effective dielectric constant
vp = c/epsilon_r # phase velocity in the medium.

### Define frequency 
pts_freq = 5000
freq = np.linspace(1e9,10e9,pts_freq)
omega = 2*np.pi*freq


# Define propagation constant
beta = omega/vp # Propagation constant. 


f0 = 5e9
lambda0 = vp/f0



def fun_open_stub(freq,L,Zs,vp = vp):
	"""Returns the Yin of open circuited stub"""
	omega = 2*np.pi*freq
	beta = omega/vp # Propagation constant. 
	#L = lambda0/4
	Zs = 50
	Yin = 1j*(1/Zs)*np.tan(beta*L)
	return Yin


def fun_photoconductive_stub(freq,L1,l,G,Zs,vp = vp):
	"""Returns the Yin of open circuited stub connected by photoconductive switch"""
	omega = 2*np.pi*freq
	beta = omega/vp # Propagation constant. 
	#L = lambda0/4
	Zs = 50
	Ys = 1/Zs
	Yin = Ys*( G*np.sin(beta*(L1+l)) + 1j*Ys*np.sin(beta*l)*np.sin(beta*L1) )/(Ys*np.sin(beta*l)*np.cos(beta*L1) - 1j*G*np.cos(beta*(L1+l))) 
	#Yin = 1j*(1/Zs)*np.tan(beta*L)
	return Yin


Z0 = 50
Zs0 = 25
ls = lambda0/10
d = 0.5*lambda0



### Define a unit cell
Cs = 0.1e-12 # Series capacitance. 

G = 0.0
Ys = 1j*Cs*omega + G

nw_gap = md.network.from_series_Z(1/Ys)

Zs0 = 10

d = 0.5*lambda0

#ls = 0.5*lambda0

#nw_tx_reso = nw_gap*mw_net.Network_Tx_line(l=ls,Z0=Zs0,gamma=1j*beta)

nw_tx_line = md.network.from_Tx_line(l=d/2,Z0=Z0,gamma=1j*beta)



nw_unit_cell = md.network.from_Tx_line(d/2,Z0,gamma=1j*beta)*nw_gap*md.network.from_Tx_line(d/2,Z0,gamma=1j*beta)


num_period = 6 # number of periods in the periodic structure
NW_filter = nw_unit_cell ** num_period # Cascade unit cell num_points times. 

#NW_filter_2 = nw_unit_cell

#for itr in range(num_period-1):
	#NW_filter_2 = NW_filter_2*nw_unit_cell



fig3 = plt.figure(3)
fig3.clf()
ax1_f3 = fig3.add_subplot(211)
ax2_f3 = fig3.add_subplot(212)

ax1_f3.plot(freq/1e9, np.abs(NW_filter.S21),label='S21')
ax1_f3.plot(freq/1e9, np.abs(NW_filter.S11),label='S11')

ax2_f3.plot(freq/1e9, 20*np.log10(np.abs(NW_filter.S21)),label='S21')
#ax2_f3.plot(freq/1e9, 20*np.log10(np.abs(NW_filter_2.S11)),label='S21')
ax2_f3.plot(freq/1e9, 20*np.log10(np.abs(NW_filter.S11)),label='S11')



ax1_f3.legend()
ax1_f3.set_xlabel('Freq (GHz)')
ax1_f3.set_ylabel('S-parameter')
ax1_f3.grid(1)

ax2_f3.legend()
ax2_f3.set_xlabel('Freq (GHz)')
ax2_f3.set_ylabel('S-parameter')
ax2_f3.set_ylim([-30,0])
ax2_f3.grid(1)

plt.show()






