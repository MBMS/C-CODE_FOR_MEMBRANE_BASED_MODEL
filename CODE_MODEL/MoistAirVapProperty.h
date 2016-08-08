// MoistAirVapProperty.h	Written by Zhiming Gao (gaoz@ornl.gov)
// Modification and extension from a code obtained from Eckhard Groll at Purdue university, used with permission.
//=================================================================================
// Unless otherwise stated, all units are MKS, with the exception of pressures (MPa) 
// mass (kg) and energy terms (kJ or KW).


#ifndef MOIST_AIRVAP_PROPERTY_INCLUDE
#define MOIST_AIRVAP_PROPERTY_INCLUDE

#include "HXM.h"
#include "SegmentProperty.h"


namespace MoistAirVapProperty
{
	double airh(double t, double rh);
	// Function airh calculates the enthalpy of moist air as a function of temperature
	// and relative humidity at atmospheric pressure.  The valid range of temperatures
	// for airh is 213.15 <= t <= 343.15 K.	[kJ/kg]

	double aircp(double t, double rh);
	// Function aircp calculates the specific heat of moist air as a function of temperature
	// and relative humidity.  The valid range of temperatures for aircp is 213.15 <= t <= 
	// 343.15 K. [kJ/kg-K]

	double airrelhum(double t, double w);
	// Function relhum calculates the relative humidity of moist air as a function of 
	// temperature and humidity ratio.  The valid range of temperatures for relhum is 213.15
	// <= t <= 343.15 K.	[-]

	double airhumrat(double t, double rh);
	// Function humrat calculates the humidity ratio of moist air as a function of 
	// temperature and relative humidity.  Data for curve fit taken from ASHRAE Fundamentals 
	// Handbook, 1997.	[kg H2O/kg dry air]

	double airmu(double t);
	// Function airmu calculates the viscosity of moist air as a function of temperature.
	// Data for curve fit taken from ASHRAE Fundamentals Handbook, 1997.	[Pa-s]

	double airk(double t);
	// Function airk calculates the thermal conductivity of moist air as a function of
	// temperature.  Data for curve fit taken from ASHRAE Fundamentals Handbook, 1997.
	// NOTE:  Use only with enthalpies given by the ASHRAE routines (e.g. airh). [kW/m-K]
	
	double airpr(double t, double rh);
	// Function airpr calculates the Prandtl number for MoistAir. [-]

	double airthsat(double h_sat);
	// Function airthsat calculates the temperature of moist air given the saturation 
	// enthalpy.	[K]

	double airthdry(double h_dry);
	// Function airthdry calculates the temperature of moist air given the dry state 
	// enthalpy.	[K]

	double airrho(double t);
	// Function airrho calculates the density of moist air as a function of temperature.
	// Data for curve fit taken from ASHRAE Fundamentals Handbook, 1997. [kg/m^3]

	double airthrh(double h, double rh);
	// Function airthrh calculates the temperature of moist air as a function of enthalpy and
	// relative humidity using bisection.	[K]

	double airthw(double h, double w);
	// Function airthw calculates the temperature of moist air as a function of enthalpy and
	// humidity ratio using bisection.	[K]

	double airrhth(double t, double h);
	// Function airrhth calculates the relative humidity of moist air given its temperature
	// and enthalpy using bisection.	[-]

	double airgetcs(double t_sat);
	// Function getcs calculates the saturated specific heat for moist air as defined in
	// Braun, et al., "Effectiveness models for cooling towers and cooling coils.", ASHRAE
	// Transactions. as a function of temperature at atmospheric pressure.	[kJ/kg-K]

	double airtdp(double t, double rh);
	// Function tdp calculates the dew point of moist air as a function of dry bulb
	// temperature and relative humidity. [K]

	double airdiff(double t);
	// Function airdiff calculates the diffusion coefficient of moisture air as a function of 
	// dry bulb temperature. [m^2/S]

	double airre_flatplate(double t, HXM* hExchanger, int flowside); 
	// Function airre_flatplate calculates the Reynolds number for moist air.  The characteristic length is 
	// demtermined based on flowside: feed side or permeate side.	[-]

	double airffr_flatplate(SegmentProperty &segpoint, HXM *hExchanger, int flowside);		
	// Function ffr_flatplate calculates the air-side friction factor for a flat-plate heat
	// exchanger surface using the correlation of Kim and Bullard.	[-]

	double airhtc_flatplate(SegmentProperty &segpoint, HXM *hExchanger, int flowside);		
	// Function htc_flatplate calculates the  air-side heat transfer coefficient for a flat-plate
	// heat exchanger surface using the correlation of Dittus-Boelter.	[kW/m^2-K]

	double airmtc_flatplate(SegmentProperty &segpoint, HXM *hExchanger, int flowside);		
	// Function mtc_flatplate calculates the  air-side mass transfer coefficient for a flat-plate
	// heat exchanger surface using the correlation of similar Dittus-Boelter.	[m/s]

	void allprop_flatplate(SegmentProperty &segpoint, HXM *hExchanger, int flowside);		
	// Function allprop_flatplate calculates all air property and htc of the  air-side for a flat-plate
	// heat exchanger surface. The reuslts are stored in SEGMENTPROPERTY struct.	[-]

	double vaph(double t);
	// Function vaph calculates the enthalpy of water vapor as a function of temperature
	// The valid range of temperatures for airh is 280K <= t <= 330 K.	[kJ/kg]

	double vapcp(double t);
	// Function vapcp calculates the specific heat of water vapor as a function of temperature
	// The valid range of temperatures for aircp is 280K <= t <= 330 K. [kJ/kg-K]

	double vapmu(double t);
	// Function airmu calculates the viscosity of water vapor as a function of temperature.
	// Data for curve fit taken from ASHRAE Fundamentals Handbook, 1997.	[Pa-s]

	double vapk(double t);
	// Function airk calculates the thermal conductivity of water vapor as a function of
	// temperature.  [kW/m-K]

	double vappr(double t);
	// Function vappr calculates the Prandtl number for water vapor.	[-]

	double vaprho(double t);
	// Function vaprho calculates the density of water vapor as a function of temperature.
	// Data for curve fit taken from ASHRAE Fundamentals Handbook, 1997.	[kg/m^3]

	double vapdiff(double t);
	// Function airdiff calculates the diffusion coefficient of water vapor as a function of 
	// dry bulb temperature. [m^2/S]

	double vapPstaturate(double t);
	// Function waterPstaturte calculates saturation pressure for the water vapor using data 
	// from the ASHRAE Fundamentals Handbook [Pa]

	double waterLatentHeat(double t);
	// Function waterLatentHeat calculates h_fg for the water that condenses out of the
	// moist air stream in the evaporator using data from the ASHRAE Fundamentals Handbook	[kJ/kg]

	double vapre_flatplate(double t, HXM* hExchanger);
	// Function airre_flatplate calculates the Reynolds number for vapor only.  The characteristic length is 
	// demtermined based on flowside: feed side or permeate side.	[-]

	double vapffr_flatplate(SegmentProperty &segpoint, HXM *hExchanger);
	// Function ffr_flatplate calculates the vapor-only permeate-side friction factor for a flat-plate heat
	// exchanger surface using the correlation of Kim and Bullard.	[-]

	double vaphtc_flatplate(SegmentProperty &segpoint, HXM *hExchanger);
	// Function htc_flatplate calculates the vapor-only permeate-side heat transfer coefficient for a flat-plate
	// heat exchanger surface using the correlation of Dittus-Boelter.	[kW/m^2-K]

	double vapmtc_flatplate(SegmentProperty &segpoint, HXM *hExchangere);
	// Function mtc_flatplate calculates the vapor-only permeate-side mass transfer coefficient for a flat-plate
	// heat exchanger surface using the correlation of similar Dittus-Boelter.	[m/s]
}

#endif
