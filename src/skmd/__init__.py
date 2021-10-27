import skmd.circuit
import skmd.data
import skmd.extract
import skmd.network
import skmd.plot
import skmd.structure 
import skmd.signals

EPSILON_0 = 8.8541878128e-12 # F/m
MU_0 = 1.25663706212e-6 # H/m
VELOCITY_OF_LIGHT = 299792458.0
OPEN = network.np.finfo(float).max
#SHORT = 1/network.np.finfo(float).max
SHORT = network.np.finfo(float).tiny
EPS = network.np.finfo(float).eps


MILLI = 1e-3
MICRO = 1e-6
NANO = 1e-9
PICO = 1e-12

GIGA = 1e9
MEGA = 1e6
KILO = 1e3
