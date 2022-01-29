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
filename4 = 'L_18pt741.csv' # S-paramters of microstripline. 


NW1 = md.data.load_HFSS_CSV(filename1)
NW2 = md.data.load_HFSS_CSV(filename2)
NW3 = md.data.load_HFSS_CSV(filename3)
NW4 = md.data.load_HFSS_CSV(filename4)



NW_gap_1 = NW4>>NW1<<NW4

Tx_parameters = md.extract.Tx_line_from_NW(NW4,l=18*md.MILLI,omega=NW4.omega)

print(Tx_parameters['L_by_l'])

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


fig2 = plt.figure('Extracted_Tx_line_parameters')
ax1_f2 = fig2.add_subplot(111)

ax1_f2.plot(NW4.freq/md.GIGA,Tx_parameters['Z0'])

plt.show()
