import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl


import skmd as md


mpl.rc('axes',labelsize=24)
mpl.rc('font',size=24)

########################################################################

filename1 = 'c_0pt25_l_18pt741.csv'
filename2 = 'c_0pt12_l_18pt741.csv'
filename3 = 'c_0pt4_l_18pt741.csv'
filename4 = 'L_18pt741.csv'


NW1 = md.data.load_HFSS_CSV(filename1)
NW2 = md.data.load_HFSS_CSV(filename2)
NW3 = md.data.load_HFSS_CSV(filename3)
NW4 = md.data.load_HFSS_CSV(filename4)



### Plot commands

fig1 = plt.figure('Loaded network')
ax1_f1 = fig1.add_subplot(111)

#ax1_f1.plot(NW1.freq/md.GIGA,md.dB_mag(NW.S11),label='$S_{11}$',linewidth=2)
#ax1_f1.plot(NW1.freq/md.GIGA,md.dB_mag(NW.S21),label='$S_{21}$',linewidth=2)
ax1_f1.plot(NW1.freq/md.GIGA,md.dB_mag(NW1.S21),label=filename1,linewidth=2)
ax1_f1.plot(NW2.freq/md.GIGA,md.dB_mag(NW2.S21),label=filename2,linewidth=2)
ax1_f1.plot(NW3.freq/md.GIGA,md.dB_mag(NW3.S21),label=filename3,linewidth=2)
ax1_f1.plot(NW4.freq/md.GIGA,md.dB_mag(NW4.S21),label=filename4,linewidth=2)
ax1_f1.legend()
plt.show()
