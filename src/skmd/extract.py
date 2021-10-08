"""Module on extraction of network parameters"""
import numpy as np


def Tx_line_Z0_gamma_from_RLCG(R,L,C,G,omega):
	"""Computes the characteristic impedance and propagation constant gamma = alpha + j beta
	Input:-
	R := Series resistance per unit length
	G := Shunt resistance per unit length
	L := Series inductance per unit length
	C := Shunt capacitance per unit length
	omega := Frequency in radians per second 
	-----------------------------------------
	Output:-
	Z0 := Characteristic impedance
	gamma := Propagation constant
	"""
	Z0 = np.sqrt((R+1j*omega*L )/(G+1j*omega*C)) 
	gamma = np.sqrt((R+1j*omega*L )*(G+1j*omega*C))
	return Z0, gamma


