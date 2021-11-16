import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#import microwave_networks as mw_net

import skmd as md


mpl.rc('axes',labelsize=24)
mpl.rc('font',size=24)

### Define constant


c = md.VELOCITY_OF_LIGHT


epsilon_r = 10 # dielectric constant or the effective dielectric constant
vp = c/epsilon_r # phase velocity in the medium.

### Define frequency 
pts_freq = 500
freq = np.linspace(0e9,10e9,pts_freq)
omega = 2*np.pi*freq


# Define propagation constant
beta = omega/vp # Propagation constant. 



#lambda0 = vp/f0


#def fun_H_omega(NW,Rd,ZL):
	#H_omega = 1/(NW.A) #+ NW.B/ZL + Rd*(NW.C + NW.D/ZL) )
	#return H_omega

def transfer_function_interconnect(NW,Rd,ZL):
	H_omega = 1/((NW.A) + NW.B/ZL + Rd*(NW.C + NW.D/ZL) )
	return H_omega


#R = 8.829e-3 * 1e6  # Ohm/meter
#G = 0
#L = 1.538e-12 * 1e6 # H/meter
#C = 0.18e-15 * 1e6 # F/meter
R = 1 # Ohm/meter
G = md.EPS
L = 0.01e-6 # H/meter
C = 100e-12 # F/meter

f0 = 4*md.GIGA

Z0_f0, gamma_f0 = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,2*np.pi*f0)

beta_f0 = gamma_f0.imag
lambda_g = 2*np.pi/beta_f0

l_inter = 0.5*lambda_g # Length of interconnect in meters

CL = 50e-15 # Load capacitanc in Farad. 

RL = 1/md.EPS
ZL = RL + 1/(1j*omega*CL + md.EPS)

Rd = 2 # Souorce resistance in Ohm. 



### Check Tx_line_Z0_gamma_from_RLCG

print('Z0 at f0 = ', Z0_f0)
print('gamma at f0 = ', gamma_f0)

Z0_omega, gamma_omega = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,omega)


#print(Z0_omega)
#print(gamma_omega)


### Define transmission line characteristic impedance and propagation constant. 


NW_interconnect = md.network.from_Tx_line(l_inter,Z0_omega,gamma_omega,omega)

H_inter = transfer_function_interconnect(NW_interconnect,Rd,ZL = ZL)



#### Along the length on interconnect. 



N_seg = 103 # Number of segments in which the interconnect is divided. 
delta_x = l_inter/N_seg

x = l_inter/N_seg*np.arange(N_seg)

#print(x)
#print(lambda_g)

def compute_load_impedance(omega,RL=RL,CL=CL):
	ZL = RL + 1/(1j*omega*CL + md.EPS)
	#ZL = np.sqrt((R+1j*omega*L)/(G+1j*omega*C))
	
	return ZL

	


def compute_H_x_omega(x,l_inter,R,L,C,G,omega):
	
	Z0_omega, gamma_omega = md.extract.Tx_line_Z0_gamma_from_RLCG(R,L,C,G,omega)


	NW_x_omega = []
	ZL_x_omega = []
	H_x_omega =[]

	for x_count in x:
		NW_count = md.network.from_Tx_line(x_count,Z0_omega,gamma_omega,omega)
		
		NW_x_omega.append(NW_count)
		
		l_minus_x = l_inter - x_count
		
		
		NW_l_minus_x = md.network.from_Tx_line(l_minus_x,Z0_omega,gamma_omega,omega)
		
		ZL_count = NW_l_minus_x.input_impedance(ZL = compute_load_impedance(omega))
		#ZL_count = NW_l_minus_x.input_impedance(ZL = 10)
		ZL_x_omega.append(ZL_count)
		
		H_x_omega.append(transfer_function_interconnect(NW_count,Rd=Rd,ZL=ZL_count))
		#H_x_omega.append(transfer_function_interconnect(NW_count,Rd=Rd,ZL=10))
		
	ZL_x_omega = np.array(ZL_x_omega)
	H_x_omega = np.array(H_x_omega)
	return H_x_omega, ZL_x_omega
		

H_x_omega , ZL_x_omega = compute_H_x_omega(x,l_inter,R,L,C,G,omega)
	
print(type(H_x_omega))
print(np.shape(H_x_omega))
	
	
#H_x = transfer_function_interconnect(NW_x,Rd, ZL = ZL_x)




#### Plot commands

#fig1 = plt.figure('Transfer_function of the interconnect')
#ax1_f1 = fig1.add_subplot(111)

#ax1_f1.plot(NW_interconnect.omega/(2*np.pi*1e9),np.abs(H_inter))
##ax1_f1.plot(NW_interconnect.omega/(2*np.pi*1e9),np.abs(NW_interconnect.S21))
#ax1_f1.set_xlabel('Freq')



fig2 =plt.figure('Properties along length')
ax1_f2 = fig2.add_subplot(211)
ax2_f2 = fig2.add_subplot(212)

index_freq = 100
print('index freq',freq[index_freq]/1e9)
ax1_f2.plot(x/lambda_g, np.abs(ZL_x_omega[:,index_freq]),label='Freq = $\\frac{\omega_1}{2\pi} = $'+str(np.round(freq[index_freq]/md.GIGA,3))+' GHz')
ax1_f2.grid()
ax1_f2.legend()

x_grid, freq_grid = np.meshgrid(omega/(2*np.pi*md.GIGA),x/lambda_g)
ax2_f2.contourf(x_grid,freq_grid, np.abs(H_x_omega))

ax1_f2.set_xlabel('$x/\lambda_g$')
ax1_f2.set_ylabel('$H(x,\omega_1)$')

ax2_f2.grid()
#fig2.tight_layout()
#ax2_f2.legend()


########################################################################
#################### Time domain analysis ##############################

"""We will have to approximate a continuous-time signal. So we divide the time 
into pts_t number of points. Larger is the 'pts_t', more accurate will be the continuous-time representation"""

pts_t = 1002 # 
T_range = 10*md.NANO

t = np.linspace(-T_range/2,T_range/2,pts_t)

vin_t = md.signals.unit_step(t-0.5*md.NANO)-md.signals.unit_step(t-2.5*md.NANO)


Vin_omega, omega_t = md.signals.Fourier_Transfrom(vin_t,t)

H_x_omega , ZL_x_omega = compute_H_x_omega(x,l_inter,R,L,C,G,omega_t)

H_omega = H_x_omega[100,:]

vout_x_t = []
Vout_omega = Vin_omega*H_omega

vout_t, t_ifft = md.signals.Inverse_Fourier_Transform(Vout_omega,omega_t)

for count_x in range(N_seg):
	Vout_omega_count = Vin_omega*H_x_omega[count_x,:]

	vout_t_count, t_ifft = md.signals.Inverse_Fourier_Transform(Vout_omega_count,omega_t)
	
	vout_x_t.append(vout_t_count)
	
	
vout_x_t = np.array(vout_x_t)


fig3 = plt.figure('Time-domain')
ax1_f3 = fig3.add_subplot(211)
ax2_f3 = fig3.add_subplot(212)

ax1_f3.plot(t/md.NANO,vin_t,linewidth=3,label='$v_{in}(t)$')
ax1_f3.plot(t_ifft/md.NANO,vout_t,linewidth=3,label='$v_{out}(t)$')
#ax1_f3.plot(t_ifft,np.real(vout_x_t[0,:]),linewidth=3,label='$v_{out}(t)$')



ax1_f3.set_xlabel('Time')
ax1_f3.grid()
ax1_f3.legend()

x_grid, t_grid = np.meshgrid(t/md.NANO,x)
#ax2_f3.contourf(np.real(vout_x_t),levels=100,cmap='hot')
ax2_f3.contourf(x_grid, t_grid, np.real(vout_x_t),levels=100,cmap='hot')
#ax2_f3.plot(omega_t/(2*np.pi*md.GIGA),np.abs(Vin_omega),linewidth=3,label='$V_{in}(\omega)$')
#ax2_f3.plot(omega_t/(2*np.pi*md.GIGA),np.abs(Vout_omega),linewidth=3,label='$V_{out}(\omega)$')
ax2_f3.set_xlabel('Time (ns)')
ax2_f3.set_ylabel('Dist')
ax2_f3.grid()
#fig3.tight_layout()
#ax2_f3.legend()

####################### Animate #############################
flag_animate = 1
flag_save = 0

if flag_animate == 1:
	fig4 = plt.figure('Animate-interconnect',figsize=(100,40))
	ax1_f4 = fig4.add_subplot(132)
	ax2_f4 = fig4.add_subplot(131)
	ax3_f4 = fig4.add_subplot(133)

	v_cmap_min = -1.5
	v_cmap_max = 1.5

	use_colormap = 'jet'

	ax2_f4.plot(t/md.NANO,vin_t,linewidth=3,label='$v_{in}(t)$')
	ax2_f4.plot(t/md.NANO,np.real(vout_x_t[0,:]),linewidth=3,label='$v(0,t)$')
	ax3_f4.plot(t/md.NANO,np.real(vout_x_t[-1,:]),linewidth=3,label='$v(l,t)$')

	md.plot.plot_colored_line(t/md.NANO,np.real(vout_x_t[0,:]),np.real(vout_x_t[0,:]),ax2_f4,vmin=v_cmap_min,vmax=v_cmap_max,use_colormap=use_colormap)

	md.plot.plot_colored_line(t/md.NANO,np.real(vout_x_t[-1,:]),np.real(vout_x_t[-1,:]),ax3_f4,vmin=v_cmap_min,vmax=v_cmap_max,use_colormap=use_colormap)

	marker_line_input, = ax2_f4.plot([],[],linewidth=2)
	marker_line_output, = ax3_f4.plot([],[],linewidth=2)

	ax2_f4.legend()
	ax3_f4.legend()

	ax2_f4.grid('on')
	ax3_f4.grid('on')

	ax2_f4.set_xlabel('Time (ns)')
	ax3_f4.set_xlabel('Time (ns)')
	ax2_f4.set_ylabel('Voltage')
	ax3_f4.set_ylabel('Voltage')

	plt.tight_layout()
	plt.pause(0.5)
	dummy = 0



	for count in range(int(pts_t/2),pts_t,5):
		y = np.zeros_like(x)


		ax1_f4.cla()
		ax1_f4.set_xlabel('Distance along the length of the interconnect',fontsize=14)
		ax1_f4.set_title('Time = %1.3f ns' %(t[count]/md.NANO))
		md.plot.plot_colored_line(x,y,np.real(vout_x_t[:,count]),ax1_f4,linewidth=10,vmin=v_cmap_min,vmax=v_cmap_max,use_colormap=use_colormap)
		marker_line_input.set_data([t[count]/md.NANO,t[count]/md.NANO],[-2,2])
		marker_line_output.set_data([t[count]/md.NANO,t[count]/md.NANO],[-2,2])
		ax1_f4.set_frame_on(False)
		ax1_f4.axes.get_yaxis().set_visible(False)
		ax1_f4.grid('on')
		if flag_save == 1:
			fname = 'interconnect-%04d.png'%dummy
			dummy = dummy+1
			print(fname)
			fig4.savefig(fname)
		plt.draw()
		plt.tight_layout()
		plt.pause(0.0001)



plt.show()
