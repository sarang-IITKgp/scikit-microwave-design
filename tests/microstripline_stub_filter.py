import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl


import skmd as md


mpl.rc('axes',labelsize=24)
mpl.rc('font',size=24)

### Define constant

c = md.VELOCITY_OF_LIGHT

### Define frequency 
pts_freq = 1000
freq = np.linspace(0.5*md.GIGA,3*md.GIGA,pts_freq)
omega = 2*np.pi*freq



#### define substrate
epsilon_r = 10.8 # dielectric constant or the effective dielectric constant
h_subs = 1.27*md.MILLI # meters. 

############### 

f0 = 1.5*md.GIGA

omega0 = md.f2omega(f0)

msl_Tx1 = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=1.1*md.MILLI,l=5*md.MILLI,text_tag='Left-line',omega=omega)
msl_Tx2 = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=1.1*md.MILLI,l=5*md.MILLI,text_tag='Right-line',omega=omega)

msl_Tx1.print_specs()

w_stub = 1*md.MILLI

lambda_g_stub = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=w_stub,text_tag='stub',omega=omega0).lambda_g



w_stub = 1*md.MILLI

L_stub = lambda_g_stub/4

msl_stub = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=w_stub,l=L_stub,text_tag='stub',omega=omega)




def define_NW_for_stub(msl_stub,ZL_stub):
	Y_stub = msl_stub.NW.input_admittance(1/ZL_stub)
	NW_stub = md.network.from_shunt_Y(Y_stub)
	return NW_stub
	

ZL_stub = md.OPEN
NW_stub = md.network.from_shunt_Y(msl_stub.NW.input_admittance(1/ZL_stub))




NW_filter = msl_Tx1.NW*NW_stub*msl_Tx2.NW



#print()

#msl_port.print_specs()
#msl_cap.print_specs()
#msl_ind.print_specs()





#msl2 = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=1.1*md.MILLI,text_tag = 'The second line')
#msl2.print_specs()


## Plot commands

fig1 = plt.figure('LPF')
ax1_f1 = fig1.add_subplot(111)
ax1_f1.plot(omega/(2*np.pi*md.GIGA),np.abs(NW_filter.S11),linewidth='3',label='$|S_{11}|$')
ax1_f1.plot(omega/(2*np.pi*md.GIGA),np.abs(NW_filter.S21),linewidth='3',label='$|S_{21}|$')
ax1_f1.grid(1)

ax1_f1.legend()

fig2 = plt.figure('LPF-mag-phase')
ax1_f2 = fig2.add_subplot(211)
ax1_cmap_f2 = fig2.add_axes([0.92, 0.1, 0.02, 0.7])
ax2_f2 = fig2.add_subplot(212)

md.plot.plot_colored_line(md.omega2f(omega)/md.GIGA,np.abs(NW_filter.S11),np.angle(NW_filter.S11)*180/np.pi,ax=ax1_f2,color_axis = ax1_cmap_f2)
md.plot.plot_colored_line(md.omega2f(omega)/md.GIGA,np.abs(NW_filter.S21),np.angle(NW_filter.S21)*180/np.pi,ax=ax2_f2)

ax1_f2.grid(1)
ax2_f2.grid(1)


fig3 = plt.figure('LPF-mag-phase-dB')
ax1_f3 = fig3.add_subplot(211)
ax1_cmap_f3 = fig3.add_axes([0.92, 0.5, 0.02, 0.3])
ax2_f3 = fig3.add_subplot(212)
ax2_cmap_f3 = fig3.add_axes([0.92, 0.1, 0.02, 0.3])

md.plot.plot_colored_line(md.omega2f(omega)/md.GIGA,md.dB_mag(NW_filter.S21),np.angle(NW_filter.S21)*180/np.pi,ax=ax1_f3,color_axis = ax1_cmap_f3)
md.plot.plot_colored_line(md.omega2f(omega)/md.GIGA,np.rad2deg(np.angle(NW_filter.S21)),md.dB_mag(NW_filter.S21),ax=ax2_f3,color_axis = ax2_cmap_f3)

ax1_f3.grid(1)
ax2_f3.grid(1)


fig4 = plt.figure('Smith-chart')


ax1_f4 = md.plot.plot_smith_chart(md.omega2f(omega)/md.GIGA,NW_filter.S21,fig4,use_colormap='inferno',linewidth=10)

snap_cursor_2 = md.plot.SnaptoCursor_polar(ax1_f4,md.omega2f(omega), NW_filter.S21)
fig4.canvas.mpl_connect('motion_notify_event', snap_cursor_2.mouse_move)




plt.show()








