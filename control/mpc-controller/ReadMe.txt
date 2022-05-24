************************************************************

			  University of Bath
		Mechanical Engineering Department
		     2022 GBDP TBRe-AI Control 


************************************************************

This file contains the the scripts and functions used to
develop the MPC controller for use on the 2022 TBRe-AI car.

To run the simulations the CasADi file call needs to be 
changed to the location the file has been saved on the users
computer. 

Controller_Track and Controller_Track_Stop_Start are the two
simulations used to demonstrate the capabilities of MPC 
around a track with the latter performing from stationary
with the necessary adjustments to the vehicle model mentioned
in the final report. The two controllers are currently tuned 
for the given vehicle parameters, tire coefficients and input
maximum values. Some of these values were not finalized at 
the time of tuning so will need to be changed dependent on the
final car.

As mentioned in the final report the MPC controller needs to 
be converted to C code and be altered for interface with
the other subsystems before use. To help with this process
a seperate script Known_Track has been started which contains
comments and adjustments to assist in this integration. It is
reccomended that the error handling code be added after 
conversion to C and integration into ROS.

************************************************************