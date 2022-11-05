"""Module containing the basic Microwave network object"""  

import matplotlib.pyplot as plt 
from matplotlib import cm
import numpy as np
import math

class Network:
	
	def __init__(self,p11,p12,p21,p22,freqpts,Nhar,parameter='abcd',Z0=50,omega=None):
		
		"""For conversion from ABCD <-> S, it is assumed that the network has equal and real characteristic impedance on both the ports."""
		
		self.Z0=Z0
		self.omega = omega
		self.freqpts = freqpts;
		self.Nhar = Nhar;
		U=np.identity(2*self.Nhar+1);
		
		self.j=2*Nhar+1;self.k=2*Nhar+1;
		
		self.A = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		self.B = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		self.C = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')	
		self.D = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		
		if parameter == 'abcd' or parameter == 'ABCD':
			
			for f in range(0,self.freqpts,1):
		
				self.A[f] = p11[f]; 
				self.B[f] = p12[f]; 	
				self.C[f] = p21[f]; 	
				self.D[f] = p22[f]; 
		

				#print('ABCD defined')
			
			self.ABCD_to_S()
			
			
	def ABCD_to_S(self,):
		
		U=np.identity(2*self.Nhar+1);
		self.S11 = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		self.S12 = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		self.S21 = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		self.S22 = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		self.S11_centre = np.zeros((self.freqpts), dtype='complex_');
		self.S21_centre = np.zeros((self.freqpts), dtype='complex_');
		self.S12_centre = np.zeros((self.freqpts), dtype='complex_');
		self.S210m2 = np.zeros((self.freqpts), dtype='complex_');
		self.S210p2 = np.zeros((self.freqpts), dtype='complex_');
		self.S210m1 = np.zeros((self.freqpts), dtype='complex_');
		self.S210p1 = np.zeros((self.freqpts), dtype='complex_');
		
		for f in range(0,self.freqpts,1):
		
			self.S11[f]=U-2*np.linalg.inv(U+(self.A[f]*self.Z0+self.B[f])@np.linalg.inv(self.C[f]*self.Z0*self.Z0+self.D[f]*self.Z0)); 
			self.S12[f]=2*U@np.linalg.inv(U+(self.A[f]*self.Z0+self.B[f])@np.linalg.inv(self.C[f]*self.Z0+self.D[f])/self.Z0)@(self.A[f]-(self.A[f]*self.Z0+self.B[f])@np.linalg.inv(self.C[f]*self.Z0+self.D[f])@self.C[f]); 
			self.S21[f]=(2*math.sqrt(self.Z0/self.Z0)*U)@np.linalg.inv(self.A[f]+self.B[f]/self.Z0+self.C[f]*self.Z0+self.D[f]*(self.Z0/self.Z0)); 
			self.S22[f]=U-2*U@np.linalg.inv(U+np.linalg.inv(self.A[f]*self.Z0+self.C[f]*self.Z0*self.Z0)@(self.B[f]+self.D[f]*self.Z0));
		#print('computed ABCD to S')
		
			self.S11_centre[f] = self.S11[f][self.Nhar][self.Nhar];
			self.S21_centre[f] = self.S21[f][self.Nhar][self.Nhar];
			self.S12_centre[f] = self.S12[f][self.Nhar][self.Nhar];
		#defined center
		
		
			self.S210m2[f] = self.S21[f][0][2];
			self.S210p2[f] = self.S21[f][4][2];
			self.S210m1[f] = self.S21[f][1][2];
			self.S210p1[f] = self.S21[f][3][2];
		
	def input_admittance(self,YL):  
		"""Returns the input admittance of a network at port-1, when Port-2 is terminated by YL"""
		"""YL should be of the same dimensions as A,B,C,D. """
		
		U=np.identity(2*self.Nhar+1);
		Yin = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		
		for f in range(0,self.freqpts,1):
			
			Yin[f] = (self.C[f]+(YL[f]@(self.D[f]*U)))@(np.linalg.inv((self.A[f]*U)+(YL[f]@self.B[f])));
			
			
		return Yin
		
		
		
	def input_impedance(self,ZL):
		"""Returns the input impedance of a network at port-1, when Port-2 is terminated by ZL"""
		"""ZL should be of the same dimensions as A,B,C,D. """
		
		U=np.identity(2*self.Nhar+1);
		Zin = np.zeros((self.freqpts,self.j, self.k), dtype='complex_');
		
		for f in range(0,self.freqpts,1):
		
			Zin[f] = ((self.A[f]@ZL[f])+self.B[f])@(np.linalg.inv((self.C[f]@ZL[f])+self.D[f]))
		
		return Zin
		
		
		
	def __mul__(self,other):
		
	
		A = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		B = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		C = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')	
		D = np.zeros((self.freqpts,self.j, self.k), dtype='complex_')
		
		for f in range(0,self.freqpts,1):
		
			A1=self.A[f]; B1=self.B[f]; C1=self.C[f]; D1=self.D[f];
		
			A2=other.A[f]; B2=other.B[f]; C2=other.C[f]; D2=other.D[f]; 
			
			A[f] = A1@A2 + B1@C2
			B[f] = A1@B2 + B1@D2
			C[f] = C1@A2 + D1@C2
			D[f] = C1@B2 + D1@D2 
			
		return Network(A,B,C,D,self.freqpts,self.Nhar,parameter='abcd')   
		
		
def Omega_mat(Freq,freq_m,Nhar,freqpts):
		
	a=(2*Nhar+1); 
	b=(2*Nhar+1);
	x=np.zeros((a,b));
	OmegaMat= np.zeros((freqpts,a,b));
	Omega = np.zeros(freqpts);
	Omega_M=2*np.pi*freq_m;  
	
	for w in range(0,freqpts,1):
		
		Omega[w] = 2*np.pi*Freq[w];
		
		for har in range(0,2*Nhar+1,1):
			
		
			x[har][har] = Omega[w]+(har-Nhar)*Omega_M;
			
		OmegaMat[w] = x;
	
	
	return OmegaMat	      
		

def admittance_of_Cap(freqpts,Nhar,OmegaMat,phi,zeta,C0):
		
	l=freqpts;a=(2*Nhar+1); b=(2*Nhar+1);U=np.identity(2*Nhar+1);
	
	YC= np.zeros((freqpts,a,b),dtype='complex_');Pmat = np.zeros((a,b),dtype='complex_');
	
	val1=(zeta/2)*np.exp(-phi*1j); 
	val2=(zeta/2)*np.exp(phi*1j);
	
	pp=0;qq=0;
	for pp in range(0,2*Nhar+1,1):
		
		if(qq>0 and qq<4):
			Pmat[pp][qq+1] = val1;
			Pmat[pp][qq-1] = val2;
		elif(qq==0):
			Pmat[pp][qq+1] = val1;
		else :
			Pmat[pp][qq-1] = val2;	
		
		qq=qq+1;
		
	for f in range(0,freqpts,1):
					
		YC[f] = 1j*C0*(OmegaMat[f]@(U + Pmat));
		
	return YC            
	
	
def from_Tx_line(L,Z0,OmegaMat,Vp,Nhar,freqpts):
		
	"""This function takes the following input paramters,
	l = length of the transmission line. Should be a scalar.
	Z0 = Characteristics impedance of the transmission line. Should be a scalar or an array of the same dimensions as beta.
	gamma = propagation constant of the transmission line. Should be the an array of the same dimensions as omega.
	
	The function returns the ABCD parameters of a transmission line section of length l."""
	
	j = 2*Nhar+1; k = 2*Nhar+1; U=np.identity(2*Nhar+1);

	A1= np.zeros((freqpts, j, k), dtype='complex_')
	B1= np.zeros((freqpts, j, k), dtype='complex_')
	C1= np.zeros((freqpts, j, k), dtype='complex_')
	D1= np.zeros((freqpts, j, k), dtype='complex_')
   
	for f in range(0,freqpts,1):
		
		beta=OmegaMat[f]/Vp; 
				
		Y0=1/Z0;
	
		A1[f]=np.cos(beta*L)*U; 
		B1[f]=1j*Z0*np.sin(beta*L)*U; 
		C1[f]=1j*Y0*np.sin(beta*L)*U; 
		D1[f]=np.cos(beta*L)*U; 
		
		
	
	return Network(A1,B1,C1,D1,freqpts,Nhar,parameter='abcd')
	
	
def from_series_Z(Z,OmegaMat,Nhar,freqpts):
	""" Returns a network object"""
	
	l=freqpts; j = 2*Nhar+1; k = 2*Nhar+1;U=np.identity(2*Nhar+1);O=np.zeros([2*Nhar+1,2*Nhar+1])
	A= np.zeros((l, j, k), dtype='complex_')
	B= np.zeros((l, j, k), dtype='complex_')
	C= np.zeros((l, j, k), dtype='complex_')
	D= np.zeros((l, j, k), dtype='complex_')
	
		
	for f in range(0,freqpts,1):
		
		A[f]=U; 
		B[f]=Z[f];
		C[f]=O;
		D[f]=U;
			
	
	return Network(A,B,C,D,freqpts,Nhar,parameter='abcd')
	 
	 
def from_shunt_Y(Y,OmegaMat,Nhar,freqpts):
	""" Returns a network object"""
	
	l=freqpts; j = 2*Nhar+1; k = 2*Nhar+1;U=np.identity(2*Nhar+1);O=np.zeros([2*Nhar+1,2*Nhar+1])
	A= np.zeros((l, j, k), dtype='complex_')
	B= np.zeros((l, j, k), dtype='complex_')
	C= np.zeros((l, j, k), dtype='complex_')
	D= np.zeros((l, j, k), dtype='complex_')
	
	for f in range(0,freqpts,1):
		
		A[f]=U; 
		B[f]=O;
		C[f]=Y[f];
		D[f]=U;
			
			
	return Network(A,B,C,D,freqpts,Nhar,parameter='abcd')
	 
	 
	 
	

		

		
		
		
		
		
			
		
