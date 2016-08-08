// Membrane.h	Written by Zhiming Gao (gaoz@ornl.gov)
//=================================================================================
// Unless otherwise stated, all units are MKS, with the exception of pressures (Pa) 
// mass (kg), and energy terms (kJ or KW).


#ifndef MEMBRANE_INCLUDE
#define MEMBRANE_INCLUDE

#include "HXM.h"
#include "SegmentProperty.h"
#include "MoistAirVapProperty.h"


namespace Membrane
{
	double knudsen_diffus(SegmentProperty &segment, HXM* hExchanger, double t_membrane);
	// Function knudsen_diffus calculates Knudsen-diffusion as a function of temperature [kg/m2-s-Pa]

	double molecular_diffus(SegmentProperty &segment, HXM* hExchanger, double t_membrane, double p_mvap, double p_mair);
	// Function molecular_diffus calculates molecular-diffusion as a function of temperature [kg/m2-s-Pa]

	double poise_flow(SegmentProperty &segment, HXM* hExchanger, double t_membrane, double p_mvap);
	// Function poise_flow calculates poiseuille-flow as a function of temperature [kg/m2-s-Pa]

	double watpermeatd(SegmentProperty &segment, HXM* hExchanger, double p_fa, double p_pf, double t_membrane, double p_mvap, double p_mair);
	// Function watpermeat calculates water vapor permeation [kg/m2-s]

	double watpermeatc(double p_fa, double p_pf);
	// Function watpermeat calculates water vapor permeation [kg/m2-s]

	double watpermeatdais(double p_fa, double p_pf, double rh_fa, double fl_fa, double t_fa, double d_fa);
	// Function watpermeat calculates water vapor permeation [kg/m2-s]

	double airpermeat(double p_fa, double p_pf, double w_fa, double w_pf);
	// Function watpermeat calculates water vapor permeation [kg/m2-s]
}

#endif
