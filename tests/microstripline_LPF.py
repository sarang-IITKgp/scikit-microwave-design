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

omega0 = 2*np.pi*1*md.GIGA

msl_port = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=1.1*md.MILLI,l=5*md.MILLI,text_tag='50_ohm_line',omega=omega)

msl_cap = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=4*md.MILLI,l=7.11*md.MILLI,text_tag='capacitive',omega=omega)
msl_ind = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=0.2*md.MILLI,l=9.81*md.MILLI,text_tag='inductive',omega=omega)


NW_filter = msl_port.NW*msl_ind.NW*msl_cap.NW*msl_ind.NW*msl_port.NW


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



plt.show()








