- [ ] [[circuit]]
- [ ] [[data]]
- [x] [[extract]]
- [x] [[network]]
- [x] [[plot]]
- [x] [[structure]]
- [x] [[signals]]

- [x] In optional variables, change to np.NaN to None. 


- [x] [[structure#Microstripline]]: Overload `*` operator to implement cascading, when NW is defined. 
	- NOT to be done in the current version. 
	- When operator overloading is done, the thumb rule is that the resultant object should be of the same data type. Cascading two MSL does not give a new MSL. So unless a new user defined 'structure' object is not created, which will have information about cascading, operator overloading the MSL or Gap objects does not makes sense. 
- [x] [[circuit]]: Overload `*` operator to implement cascading. 
	- Not to be done in the current version. 
	- Reason: Same as for operator overloading of objects in Structure module. 
	
-----------

## Demo

#### Structure
- [x] Microstrip line design.
- [x] Microstrip line with gap.
- [x] Filter design: Stub filter.
- [x] filter design: Band-pass filter with gap.
- [x] filter desgin: periodic structure. 
- [x] filter design: periodic structure with defect.


#### Circuit and extract

- [ ] Extract transmission line parameters
	- data, extract
- [ ] De-embed msl and extract pi-model of gap. 
	- data, extract, circuit
- [ ] Impedance matching
	- circuit#transmissionline, plot


#### Network

- [x] Interconnect: time-domain visualization
	- plot, smith chart, signals




#### Data

- [ ] Test loading and plotting of data. 


