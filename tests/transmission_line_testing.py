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
pts_freq = 1000
freq = np.linspace(1e9,3e9,pts_freq)
omega = 2*np.pi*freq


# Define propagation constant
beta = omega/vp # Propagation constant. 


f0 = 3e9
lambda0 = vp/f0


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
#R = 0.1  # Ohm/meter
#G = 0.01
#L = 0.01e-6 # H/meter
#C = 100e-12 # F/meter



#Z0_f0, gamma_f0 = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,2*np.pi*f0)

#beta_f0 = gamma_f0.imag
#lambda_g = 2*np.pi/beta_f0

l_tx = lambda0*0.45 # Length of interconnect in meters

CL = 50e-15 # Load capacitanc in Farad. 
Rd = 0 # Souorce resistance in Ohm. 


### Check Tx_line_Z0_gamma_from_RLCG



#print('Z0 = ', Z0_f0)
#print('gamma = ', gamma_f0)

#Z0_omega, gamma_omega = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,omega)

gamma =  1j*beta
Z0=50
nw_tx_line = md.network.from_Tx_line(l_tx,Z0,gamma,omega)

#Z_in = nw_tx_line.input_impedance(ZL=)
print(md.OPEN)
print(md.SHORT)

print(nw_tx_line)

Zin = nw_tx_line.input_impedance(md.OPEN)
Yin = nw_tx_line.input_admittance(md.SHORT)
ZL = 30
ref = nw_tx_line.input_reflection(ZL=30)
#Yin = nw_tx_line.C/nw_tx_line.A
#Zin_open_load = nw_tx_line.B/nw_tx_line.D


def fun_tx_along_length(l_tx,Z0,ZL,gamma,omega,N_seg=100):
	
	delta_x = l_tx/N_seg
	tx_delta = md.network.from_Tx_line(delta_x,Z0,gamma,omega)
	tx_x = []
	Z_x = []
	ref_x = []
	#tx_count = tx_delta
	for count in range(N_seg):
		tx_count = tx_delta**(count+1)
		#tx_count = tx_delta*tx_count
		tx_x.append(tx_count)
		print(count)
		Z_x.append(tx_count.input_impedance(ZL))
		print(np.shape(Z_x))
		ref_x.append(tx_count.input_reflection(ZL,Z0=Z0))
	#tx_x = np.array(tx_x)
	#Z_x = np.array(Z_x)
	#ref_x = np.array(ref_x)
	#print(type(Z_x))
	#print(np.shape(Z_x))
	#print(Z_x)
	#print(type(tx_x))
	#print(np.shape(tx_x))
	#print(tx_x)
	return tx_x, np.array(Z_x), np.array(ref_x)
	
	
	

N_seg = 50
omega0 = 2*np.pi*f0
beta0 = omega0/vp
gamma0 = 1j*beta0

l_var = l_tx/N_seg*np.arange(N_seg)/lambda0

tx_x, Z_x, ref_x = fun_tx_along_length(l_tx,Z0,ZL,gamma0,omega0,N_seg=N_seg)

fig1 = plt.figure()
ax1_f1 = fig1.add_subplot(111)
#ax1_f1.plot(nw_tx_line.freq,np.real(Zin),label='Re{Zin}')
#ax1_f1.plot(nw_tx_line.freq,np.imag(Zin),label='Im{Zin}')
#ax1_f1.plot(nw_tx_line.freq,np.real(Yin),label='Re{Yin}')
#ax1_f1.plot(nw_tx_line.freq,np.imag(1/Yin),label='1/Im{Yin}')
ax1_f1.plot(nw_tx_line.freq,np.real(ref),label='$Re{\Gamma}$')
ax1_f1.plot(nw_tx_line.freq,np.imag(ref),label='$Im{\Gamma}$')
#ax1_f1.plot(nw_tx_line.freq,np.abs(ref),label='$|\Gamma|$')
ax1_f1.legend()
#ax1_f1.plot(nw_tx_line.freq,np.abs(nw_tx_line.S11))
ax1_f1.set_ylim([-5,5])
ax1_f1.grid(1)


fig2 = plt.figure('Smith-chart')

#ax1_f2 = fig2.add_subplot(111)


ax1_f2 = md.plot.plot_smith_chart(l_var,ref_x,fig2,use_colormap='inferno',linewidth=10)
#ax1_f2 = md.plot.plot_smith_chart(freq,ref,fig2,use_colormap='inferno',linewidth=10)
#ax1_f2 = plot_smith_chart(freq,NW_filter.S21,fig2,linewidth=7,use_colormap='jet')


snap_cursor_2 = md.plot.SnaptoCursor_polar(ax1_f2,l_var, ref_x)
#snap_cursor_2 = md.plot.SnaptoCursor_polar(ax1_f2,nw_tx_line.freq/1e9, ref)
fig2.canvas.mpl_connect('motion_notify_event', snap_cursor_2.mouse_move)


#hand_dynamic = md.plot.Dynamic_data_points(fig2,ax1_f2,nw_tx_line.freq,ref)


fig3 = plt.figure('Along_the_line')
ax1_f3 = fig3.add_subplot(211)
ax2_f3 = fig3.add_subplot(212)

ax1_f3.plot(l_var, np.real(Z_x))
ax2_f3.plot(l_var, np.real(ref_x))
ax2_f3.plot(l_var, np.imag(ref_x))



plt.show()
