//Component.h written by Zhiming Gao (gaoz@ornl.gov)


#ifndef Heat_Exchanger_Component
#define Heat_Exchanger_Component

#include <cmath>
#include "HXM.h"
#include "SegmentProperty.h"
#include "MoistAirVapProperty.h"
#include "Membrane.h"
#include "SegmentSolve.h"
#include "Cmatrix.h"


// General heat exchanger component simulation based segment
class Component
{
	public:
		double length;				// L (m)
		double depth;				// D (m)
		double thick_fa;			// Hfa (m)
		double thick_pf;			// Hpf (m)
		double thick_membrane;		// Hm (m)
		double mf_fa;				// total feed airflow (kg/s)
		double mfvap_fa;			// inlet total vapr flow of feed airflow (kg/s)
		double mf_pf;				// total permeate flow (kg/s)
		double mfvap_pf;			// inlet total vapr flow of permeate flow (kg/s)
		double rp_membrane;			// radius of membrane pore, (m)
		double ps_membrane;			// membrane porosity, (dimensionless)
		double ts_membrane;			// membrane tortuosity (dimensionless)
		double kt_membrane;			// membrane thermoconductivity (kW/m-K)
		double temp_fa;				// feed air temperature (K)
		double vap_fa;				// feed air vaporflow (kg/s)
		double humid_fa;			// feed air humidity (kg/kg)
		double press_fa;			// feed air pressure (Pa)
		double flow_fa;				// feed air  mass rate (kg/s)
		double temp_pf;				// permeate flow  temperature (K)
		double vap_pf;				// permeate flow  vaporflow (kg/s)
		double humid_pf;			// permeate flow  humidity (kg/kg)
		double press_pf;			// permeate flow  pressure (Pa)
		double flow_pf;				// permeate flow  mass rate (kg/s)

		double tfa_out;				// component-out temperature at feed side (K)
		double vapfa_out;			// component-out vapor at feed side (kg/s)
		double flfa_out;			// component-out mass rate at feed side (Pa)
		//double pfa_out;			// component-out mass rate at feed side (kg/s)
		double tpf_out; 			// component-out temperature at permeate side (K)
		double vappf_out;			// component-out vapor at permeate side (kg/s)
		double flpf_out; 			// component-out mass rate at permeate side (kg/s)
		//double ppf_out; 			// component-out mass rate at permeate side (Pa)

		double htc_correction;		// htc correction
		double ftc_correction;		// ftc_correction
		double memsurf_correction;	// membrane surface correction

		int num_lsegment;			// dimensionless
		int num_dsegment;			// dimensionless
		int num_layer;				// parallel membrane layers
		int flowtype;				// 1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or other: crossflow
		int proppf_airvap;			// property for permeate flow : 0: vapor only; 1: air

		Component(double lvalue, double dvalue, double tfavalue, double tpfvalue, double tmvalue, double mffa, double mfpf, 
		  double rpmvalue, double psmvalue, double tsmvalue, double ktmvalue, 
		  double tfa_in, double vapfa_in, double pfa_in, double tpf_in, double vappf_in, double ppf_in, 
		  int prop_pfav, int prop_memmodel, int lnvalue, int dnvalue, int laynvalue, int ftype, double htc_correct, double ftc_correct,
		  double * *map_temp_fa, double * *map_temp_fam, double * *map_vap_fa, double * *map_humid_fa, double * *map_flow_fa, double * *map_pvap_fa,
		  double * *map_temp_pf, double * *map_temp_pfm, double * *map_vap_pf, double * *map_humid_pf, double * *map_flow_pf, double * *map_pvap_pf, double * *map_permeat_vap,
		  double enable_mtimpact, double memsurf_correct);


	private:
		void SgtPropUpdate(SegmentProperty &segpoint, HXM *hExchanger);	
};

#endif