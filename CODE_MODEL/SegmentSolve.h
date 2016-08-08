//SegmentSolve.h written by Zhiming Gao (gaoz@ornl.gov)


#ifndef SegmentSolve_H
#define SegmentSolve_H

#include <cmath>
#include "MoistAirVapProperty.h"

using namespace std;


class SegmentSolve
{
	public:
		double mflwseg_fa;		//mfl_fa: Kg/s
		double fa_rho;			//rho_fa: kg/m3
		double fa_cp;			//cp_fa: kJ/kg-k
		double lseg_fa;			//dw_fa: m
		double dseg_fa;			//dx_fa: m
		double fa_htc;			//htc_fa: kW/m2-K
		double fa_mtc;			//mtc_fa: m/s
		double fa_tin;			//T_fa_in: K
		double fa_vapin;		//vap_fa_in: Kg-vap/s
		double mflwseg_pf;		//mfl_fa: Kg/s-
		double pf_rho;			//ro_pf: kg/m3
		double pf_cp;			//cp_pf: kJ/kg-k
		double lseg_pf;			//dw_pf: m
		double dseg_pf;			//dx_pf: m
		double pf_htc;			//htc_pf: kW/m2-K
		double pf_mtc;			//mtc_pf: m/s
		double pf_tin;			//T_pf_in: K
		double pf_vapin;		//vap_pf_in: Kg-vap/s
		double delta_mem;		//thick_membrane: m
		double kt_mem;			//membrane thermal conductivity: kW/m-K
		double mfl_mem;			//membrane permeate flow rate: kg/m2-s
		double enable_mtm;		//active mt impact on ht: 0-not active; 1: active
		double enable_msf;		//active membrane deflection impact: 0-not active; 1: active
		int flowtype;			//1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or other: crossflow
		int num_membrane;		//membrane layer

		SegmentSolve();
		SegmentSolve(double mflwseg_fa_value, double fa_rho_value, double fa_cp_value, double lseg_fa_value, double dseg_fa_value, 
					 double fa_htc_value, double fa_mtc_value, double fa_tin_value, double fa_vapin_value,
					 double mflwseg_pf_value, double pf_rho_value, double pf_cp_value, double lseg_pf_value, double dseg_pf_value, 
					 double pf_htc_value, double pf_mtc_value, double pf_tin_value, double pf_vapin_value,
					 double delta_mem_value, double kt_mem_value, double mfl_mem_value, int ftype, int mem_laynvalue, 
					 double enable_mtimpact, double memsurf_correct); //initilization


		void SegmentSolht();	//solve heat transfer equation
		void SegmentSolmt();	//solve mass transfer equation

		double get_fa_temp();	//get fa_temp
		double get_fam_temp();	//get fam_temp
		double get_pfm_temp();	//get pfm_temp
		double get_pf_temp();	//get pf_temp
		double get_fa_vap();	//get fa_w
		double get_fam_vap();	//get fam_w
		double get_pfm_vap();	//get pfm_w
		double get_pf_vap();	//get pf_w
		double ask_fa_coeff();	//get fa_coeff
		double ask_pf_coeff();	//get pf_coeff
		void get_fa_coeff();	//get fa_coeff
		void get_pf_coeff();	//get pf_coeff


	private:
		double a[4][4];
		double b[4];
		double fa_temp_value, fam_temp_value, pfm_temp_value, pf_temp_value;
		double fa_vap_value, fam_vap_value, pfm_vap_value, pf_vap_value;
		double fa_coeff, pf_coeff;

		void jfc(double a[4][4], double b[4], int m, double eps);
};

#endif