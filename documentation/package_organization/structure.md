This module defines objects of the following class. 


## Class
 ###### Microstripline
 
					<----W---->
					============	^ 
									|
							er 		h
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
	
###### MSL_gap:
	Creates an object for a gap in microstrip line.
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