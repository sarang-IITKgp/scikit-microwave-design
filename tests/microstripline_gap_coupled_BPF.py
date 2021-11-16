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

w_line = 1.1*md.MILLI

msl_Tx1 = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=w_line,l=5*md.MILLI,text_tag='Left-line',omega=omega)


w_res = 1.1*md.MILLI

lambda_g_res = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=w_res,text_tag='Resonant-line-for-lambda_g',omega=omega0).lambda_g

L_res = lambda_g_res/2


msl_res = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=w_res,l=L_res,text_tag='Resonant-line',omega=omega)

d_gap = 0.125*md.MILLI

gap_in_msl = md.structure.MSL_gap(d_gap=d_gap,w=w_res, er=epsilon_r,h=h_subs,omega=omega)


#NW_filter = msl_Tx1.NW *gap_in_msl.NW*msl_res.NW*gap_in_msl.NW* msl_Tx1.NW
NW_filter = msl_Tx1.NW *gap_in_msl.NW*(msl_res.NW*gap_in_msl.NW)**3 * msl_Tx1.NW







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
md.plot.plot_colored_line(md.omega2f(omega)/md.GIGA,np.rad2deg(np.angle(NW_filter.S21)),md.dB_mag(NW_filter.S21),ax=ax2_f3,color_axis = ax2_cmap_f3,vmin=-20,vmax=0)


ax1_f3.set_ylim([-15,0])


ax1_f3.grid(1)
ax2_f3.grid(1)



plt.show()








