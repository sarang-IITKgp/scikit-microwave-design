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
freq = np.linspace(1e9,3e9,pts_freq)
omega = 2*np.pi*freq



#### define substrate
epsilon_r = 10.8 # dielectric constant or the effective dielectric constant
h_subs = 1.27*md.MILLI # meters. 

############### 

omega0 = 2*np.pi*1*md.GIGA

msl1 = md.structure.Microstripline(er=epsilon_r,h=h_subs,Z0=93,text_tag='Line-1')
msl1.print_specs()


msl2 = md.structure.Microstripline(er=epsilon_r,h=h_subs,w=1.1*md.MILLI,text_tag = 'The second line')
msl2.print_specs()

msl3 = md.structure.Microstripline(er=epsilon_r,h=h_subs,Z0=100,text_tag = 'The second line-with-frequency',omega=omega)
msl3.print_specs()













