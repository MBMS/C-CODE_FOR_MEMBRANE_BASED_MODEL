//HXM.h written by Zhiming Gao (gaoz@ornl.gov)


#ifndef Heat_Exchanger_Membrane_based
#define  Heat_Exchanger_Membrane_based

#include <cmath>

using namespace std;


// General heat exchanger properties
class HXM
{
	public:
		double length;				// L (m)
		double depth;				// D (m)
		double thick_fa;			// Hfa (m)
		double thick_pf;			// Hpf (m)
		double thick_membrane;		// Hm (m)
		double mf_fa;				// total feed airflow (kg/s)
		double mfvap_fa;			// total vapor flow in feed airflow (kg/s)
		double mf_pf;				// total permeate flow (kg/s)
		double mfvap_pf;			// total vapor flow in permeate flow (kg/s)
		double memsur_corr;			// membrane surface correction, dimensionless
		int num_lsegment;			// dimensionless
		int num_dsegment;			// dimensionless
		int num_layer;				// parallel membrane layers
		int flowtype;				// 1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or other: crossflow

		HXM();
		HXM(double lvalue, double dvalue, double tfavalue, double tpfvalue, double tmvalue, double mffa, double mfvapfa, double mfpf, double mfvappf, double memsurfr, int lnvalue, int dnvalue, int laynvalue, int ftype);

		double aseg_ht();			//  heat transfer area per segment(m2)
		double aseg_crossfa();		//  cross area of feed flow per segment (m2)
		double aseg_crosspf();		//  cross area of permeate flow per segment (m2)
		double get_lseg_fa();		//  length of feed flow per segment (m)
		double get_dseg_fa();		//  depth of feed flow per segment (m)
		double get_mflwseg_fa();	//  feed mass flow per segment (kg/s-per segment)
		double get_mflwvapseg_fa();	//  vapor flow of feed mass flow per segment (kg/s-per segment)
		double get_lcseg_fa();		//	characteristic lenght of feed airflow segment(m)
		double get_lseg_pf();		//  length of permeate flow per segment (m)
		double get_dseg_pf();		//  depth of permeate flow per segment (m)
		double get_mflwseg_pf();	//  permeate mass flow per segment (kg/s-per segment)
		double get_mflwvapseg_pf();	//  vapor flow of permeate mass flow per segment (kg/s-per segment)
		double get_lcseg_pf();		//	characteristic lenght of permeate flow segment(m)

	private:
		double lseg_fa;				// segment lenght along feed airflow (m)
		double dseg_fa;				// segment depth along feed airflow (m)
		double mflwseg_fa;			// feed airflow per segment (kg/s-per segment per layer)
		double mflwvapseg_fa;		// vapor flow of feed airflow per segment (kg/s-per segment per layer)
		double lseg_pf; 			// segment lenght  along permeate flow (m)
		double dseg_pf;				// segment depth along permeate flow (m)
		double mflwseg_pf; 			// permeateflow per segment (kg/s-per segment per layer)
		double mflwvapseg_pf; 		// vapor flow of permeateflow per segment (kg/s-per segment per layer)
		int numchannel_fa; 			// feed flow channel number along device height
		int numchannel_pf;	 		// permeate flow channel number along device height
};

#endif