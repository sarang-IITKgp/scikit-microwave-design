import skmd.circuit
import skmd.data
import skmd.extract
import skmd.network
import skmd.plot
import skmd.structure 

EPSILON_0 = 8.8541878128e-12 # F/m
MU_0 = 1.25663706212e-6 # H/m
VELOCITY_OF_LIGHT = 299792458.0
OPEN = network.np.finfo(float).max
SHORT = 1/network.np.finfo(float).max

