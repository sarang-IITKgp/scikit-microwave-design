"""Module containing physical structures"""

import numpy as np

class Microstripline:
	"""
			<----W---->
			============^ 
						|
				er		h
						|
	=============================

	Microstrip line 

	Reference: Hong & Lancaster.
	
	Defines an object of the type Microstrip line. The object can be defined in either of the following two ways. 
	
	1) Defining the characteristics impedance Z0 of the line. In this  method, the required width of the msl will be computed using synthsis equations. 
	
	2) Defining the width W of the line. In this method, the resultant characteristics impedance is computed using the analysis equations. 
	
	If both W and Z0 are given during the definition, then priority is given to W over Z0 for defining the msl. Either one of these must be defined. 
	
	Attributes:
	- Substrate dielectric constant: er (Given by user)
	
	- Substrate thickness: h (Given by user. In units of meters. Must.)
	
	- Thickness of the msl metal conductor: t (Given by the user. In units of meters. OPTIONAL. Not implemented yet.)
	
	- Effective substrate dielectric constant: er_eff (Computed)
	
	- Width of the msl: W (In units of meters. Computed from Z0, if not defined by the user.)
	
	- Characteristics impedance: Z0 (In units of Ohm. Computed if W is defined. If W is not defined, takes the value given by the user. )
	
	- If both W and Z0 are not defined: Gives an error. 
	
	- Length of the msl: l (In units of meters. Preferably should be given by the user. Takes value=1 meter, if not defined by the user.)
	
	
	- Frequency: omega. (In units of rad/sec. Must if an object of Network  class has to be computed as an attribute. Optional is only using it for synthesis and analysis.)
	
	- An object of class Network: NW (Object type: Network. Computed and defined only when omega is defined.)
	
	- Guided wavelength: lambda_g (In units of meters. Computed when omega is defined.)
	
	- Phase velocity: vp (In units of meters/sec. Computed)
	
	- Phase constant: beta (In units of rad/m. Computed when omega is defined.)
	"""
	
	def __init__(self, er,h,l=1, w=np.NaN , Z0 = np.NaN, omega=np.nan, t = np.NaN):
		self.eta = 120*np.pi # Free space impedance. In units of Ohm. 
		self.c = 3e8 # velocity of light in units of m/s.
		self.er = er
		self.h = h
		self.l = l
		#self.w = w

		self.omega = omega
		#self.t = t # Effect of conductor thickness not implemented yet. 
		
		print("============ \n Microstrip line defined")
		
		if not np.isnan(w):
			self.w = w
			print('MSL defined with width.')
			self.fun_msl_w_to_Z0()
			
		elif not np.isnan(Z0):
			print('MSL defined with Z0')
			self.Z0 = Z0
			self.fun_msl_Z0_to_w()
		
		else:
			#print("Error. Neither of W or Z0 defined.")
			sys.exit("Error in defining microstrip line ...!! Neither of W or Z0 defined. Abort...!!")
			
		self.vp = self.c/np.sqrt(self.er_eff)
		
		if not np.all(np.isnan(self.omega)):
			self.beta = self.omega/self.vp
			self.lambda_g = 2*np.pi/self.beta
			print('Frequency defined.')
			self.NW = Network_Tx_line(l=self.l,Z0=self.Z0,gamma=1j*self.beta)
		
		
	def fun_add_frequency(self,omega):
		self.omega = omega
		self.beta = self.omega/self.vp
		self.lambda_g = 2*np.pi/self.beta
		print('Frequency added (override old values)')
		self.NW = Network_Tx_line(l=self.l,Z0=self.Z0,gamma=1j*self.beta)
		
		
	def fun_msl_w_to_Z0(self,):#er,h,w):
		""" Computes Z0 and er_eff for given microstrip line parameters
		Ref: Page-49, Section-4.1.4
		"""
		u = self.w/self.h
		
		a = 1 + 1/49 * np.log((u**4 + (u/52)**2 ) /(u**4 + 0.432 ) ) \
			+ 1/18.7 * np.log(1+(u/18.1)**3 )
		
		b = 0.564*((self.er-0.9)/(self.er +3) )**0.053
		
		self.er_eff = (self.er + 1)/2 + 0.5*(self.er - 1)*(1 + 10/u)**(-a*b) # Eq 4.4
		
		F = 6 + (2*np.pi -6)*np.exp(-(30.666/u)**0.7528)
		
		self.Z0 = 120*np.pi/(2*np.pi*np.sqrt(self.er_eff)) * np.log(F/u + np.sqrt(1+(2/u)**2)) # Eq 4.5
		#return er_eff, Z0 


	def fun_msl_Z0_to_w(self,):#er,h,Z0):
		""" Computes the required with of a microstrip line for desired 
		characteristics impedance."""
		
		A = self.Z0/60 * ( (self.er+1)/2)**0.5 + (self.er-1)/(self.er+1)*(0.23 + 0.11/self.er) # Eq 4.10 '+ 1'.
		
		u = 8 * np.exp(A)/ (np.exp(2*A) - 2) # Eq. 4.10
		
		#if W_by_h < 2:
			#W = W_by_h*h
		if  u >2:
			B = 60*np.pi**2/(self.Z0*np.sqrt(self.er))
			u = 2/np.pi * ( (B-1) - np.log(2*B - 1 )  \
				+ (self.er-1)/(2*self.er)* (np.log(B-1) + 0.39 - 0.61/self.er))  # Eq. 4.11
				
		self.w = u * self.h
		
		a = 1 + 1/49 * np.log((u**4 + (u/52)**2 ) /(u**4 + 0.432 ) ) \
			+ 1/18.7 * np.log(1+(u/18.1)**3 )
		
		b = 0.564*((self.er-0.9)/(self.er +3) )**0.053
		
		self.er_eff = (self.er + 1)/2 + 0.5*(self.er - 1)*(1 + 10/u)**(-a*b) # Eq 4.4
		
		#return self.er_eff, self.W

	

#End_Microstripline_Class:



class MSL_gap:
	"""Creates an object for a gap in microstrip line.
	Computes the approximate capacitances for PI-equivalent of a MSL-gap.
					   ____
	 o----------------|_Cg_|---------------------o
			   _|__         _|__
			  |_Cp_|       |_Cp_|
				|  			 |
	o-------------------------------------------o
	Attributes include:
	- Gap length : d . In units of meters. Must be given by the user.
	- Width of the MSL: w. In units of meters. Must be given by the user.
	- Substrate thickness: h. In units of meters. Must be given by the user. 
	- Dielectric constant: er. Must be given by the user. 
	- Frequency: omega. In units of radians/second. Optional. 
	
	- Series gap capacitance: Cg. In units of Farads. Computed. 
	- Shunt capacitances: Cp. In units of Farads. Computed. 
	- Microwave network: NW. Object of the Network class. Computed only when omega is defined. 
	
	Based on Section 4.3.1.3 of the Hong & Lancaster book. 
	- Point to cross-check. Whether log in the book means log10 or natural log. In the current implementation, log10 is assumed. 
	"""
	
	def __init__(self, d_gap, w, er, h, omega=np.nan):
		
		print('Gap Class Warning: capacitance values may be in accurate')
		self.d = d_gap
		self.w = w
		self.er = er
		self.h = h
		self.omega = omega
		
		w_by_h = self.w/self.h
		
		mo = w_by_h *(0.619*np.log10(w_by_h) - 0.3853)
		ko = 4.26 - 1.453*np.log10(w_by_h)
		
		if self.d/self.w <0.3:
			me = 0.8675
			ke = 2.043*w_by_h**0.12
		else :
			me = 1.565/(w_by_h**0.16) - 1
			ke = 1.97 - 0.03/w_by_h
			
		Co = self.w*(self.er/9.6)**0.8 * (self.d/self.w)**mo * np.exp(ko) *1e-12
		Ce  = 12*self.w*(self.er/9.6)**0.9 * (self.d/self.w)**me *np.exp(ke) *1e-12
		
		self.Cp = 0.5*Ce
		self.Cg = 0.5*Co - 0.25*Ce
		
		if not np.all(np.isnan(self.omega)):
			
			self.NW = Network_PI_Y(1j*self.omega*self.Cp,1j*self.omega*self.Cp,1j*self.omega*self.Cg)
		
