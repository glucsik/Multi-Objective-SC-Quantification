# Multi-Objective-SC-Quantification
Preprocessing side of SC Quantification

Post_Process file contains an algorithm to be used in the process of
quantifying the Specific Conductance of a Stream. The Eckhardt filtering method is used to
distiguish the contributions from Baseflow (bf) and Qucikflow (qf) to total flow. 

Post_Process_v5 provides the physically based upper and lower bounds of the SC of bf and qf based on the Eckhardt's bf quantification.
These bounds are to be used as contraints for a multi objective optimization process in which the parameters, x1,x2,..., are solved for. 

The Tule_SC_RB_O_UnT_1.csv contains the Streams total flow, total SC, Temp, and PH from the USGS monitoring site at the North fork of the Tule River. 

Only total SC and Streamflow are needed to determine the SC contributions of bf and qf. 
PH and Tempurature can be included as additional means by which missing data can be interpolated.



