import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

#import microwave_networks as mw_net
import skmd as md


mpl.rc('axes',labelsize=24)
mpl.rc('font',size=24)

### Define constant
#c = 3e8 # Velocity of light in meters/second
c = md.VELOCITY_OF_LIGHT


epsilon_r = 10 # dielectric constant or the effective dielectric constant
vp = c/epsilon_r # phase velocity in the medium.

### Define frequency 
pts_freq = 5000
freq = np.linspace(1e9,100e9,pts_freq)
omega = 2*np.pi*freq


# Define propagation constant
beta = omega/vp # Propagation constant. 



#lambda0 = vp/f0


def fun_H_omega(NW,Rd,ZL):
	H_omega = 1/(NW.A) #+ NW.B/ZL + Rd*(NW.C + NW.D/ZL) )
	return H_omega


Z0 = 50
Zs0 = 25
ls = lambda0/10
d = 1*lambda0



### Define a unit cell


#R = 8.829e-3 * 1e6  # Ohm/meter
#G = 0
#L = 1.538e-12 * 1e6 # H/meter
#C = 0.18e-15 * 1e6 # F/meter
R = 0.1  # Ohm/meter
G = 0.01
L = 0.01e-6 # H/meter
C = 100e-12 # F/meter

f0 = 2e9

Z0_f0, gamma_f0 = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,2*np.pi*f0)

beta_f0 = gamma_f0.imag
lambda_g = 2*np.pi/beta_f0

l_tx = lambda_g # Length of interconnect in meters

CL = 50e-15 # Load capacitanc in Farad. 
Rd = 0 # Souorce resistance in Ohm. 


### Check Tx_line_Z0_gamma_from_RLCG



print('Z0 = ', Z0_f0)
print('gamma = ', gamma_f0)

Z0_omega, gamma_omega = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,omega)

print(Z0_omega)

#nw_tx_line = md.network.from_Tx_line(l=l_inter,Z0=Z0_inter,gamma=gamma_inter)


#ZL = 1j*omega*CL

#H_omega_interconnect = fun_H_omega(nw_tx_line,Rd,Z0_inter)

##print(size)



##nw_unit_cell = md.network.from_Tx_line(d/2,Z0,gamma=1j*beta)*nw_gap*md.network.from_Tx_line(d/2,Z0,gamma=1j*beta)


##num_period = 6 # number of periods in the periodic structure
##NW_filter = nw_unit_cell ** num_period # Cascade unit cell num_points times. 

##NW_filter_2 = nw_unit_cell

##for itr in range(num_period-1):
	##NW_filter_2 = NW_filter_2*nw_unit_cell



#fig3 = plt.figure(3)
#fig3.clf()
#ax1_f3 = fig3.add_subplot(211)
#ax2_f3 = fig3.add_subplot(212)

##ax1_f3.plot(freq/1e9, np.angle(H_omega_interconnect),label='|H(\omega)|')
#ax1_f3.plot(freq/1e9, np.abs(H_omega_interconnect),label='|H(\omega)|')
##ax1_f3.plot(freq/1e9, np.abs(NW_filter.S11),label='S11')

#ax2_f3.plot(freq/1e9, np.abs(nw_tx_line.S21),label='S21')
###ax2_f3.plot(freq/1e9, 20*np.log10(np.abs(NW_filter_2.S11)),label='S21')
##ax2_f3.plot(freq/1e9, 20*np.log10(np.abs(NW_filter.S11)),label='S11')



##ax1_f3.legend()
##ax1_f3.set_xlabel('Freq (GHz)')
##ax1_f3.set_ylabel('S-parameter')
##ax1_f3.grid(1)

##ax2_f3.legend()
##ax2_f3.set_xlabel('Freq (GHz)')
##ax2_f3.set_ylabel('S-parameter')
##ax2_f3.set_ylim([-30,0])
##ax2_f3.grid(1)

#plt.show()






