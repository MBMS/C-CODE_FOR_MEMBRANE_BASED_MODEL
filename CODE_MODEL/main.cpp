// main.cpp	written by Zhiming Gao (gaoz@ornl.gov) on 10/26/2015


// membraned-based dehumidifier system analysis program
// This is a shell main program for DAIS, software for modeling membraned-based dehumidifier systems.

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <new>
#include <typeinfo>
#include <exception>
#include <stdexcept>
#include <ios>

#include "HXM.h"
#include "SegmentProperty.h"
#include "MoistAirVapProperty.h"
#include "Membrane.h"
#include "SegmentSolve.h"
#include "Component.h"
#include "Cmatrix.h"
#include "iofile.h"

using namespace std;
using namespace Constant_Coefficients;
using namespace MoistAirVapProperty;
using namespace Membrane;
using namespace DeviceIO;


int main()
{
	try
	{
		//Read input data
		inputfile();

		//Geometry sample and allocation
		//======================================================
		//             fa_in|pf_in           : parallel uncounterflow
		//			 ______|____|____
		//	--> Hpf	|      |    v   |\ 
		//		Hm	|======v========| \
		//		Hfa	|_______________|\ \
		//			\      L        \D\|
		//			 \           |   \ |
		//			  \__________|____\|
		//				   |     V pf_out
		//				   |
		//				   V fa_out		
		//======================================================
		//                 fa_in^pf_out      : parallel counterflow
		//			 ______|____|____
		//	--> Hpf	|      |    |   |\ 
		//		Hm	|======v========| \
		//		Hfa	|_______________|\ \
		//			\      L        \D\|
		//			 \           ^   \ |
		//			  \__________|____\|
		//				   |     | pf_in
		//				   |
		//				   V fa_out
		//======================================================
		//             fa_int			     : crossflow
		//			 ______|________
		//	-->	Hpf	|      |        |\ 
		//		Hm	|======v========| \
		//		Hfa	|_______________|\-\-> pf_out
		//			\      L        \D\|
		//			 \               \ |
		//			  \_______________\|
		//					   |
		//					   |
		//					   V fa_out
		//======================================================	
		double l_dev = length_dev;					// length of device, m
		double d_dev = deep_dev;					// deep of device, m
		double hfa_dev = height_feedchanel;			// height of feedair chanel, m
		double hpf_dev = height_permeatechanel;		// height of permeateflow chanel, m
		double tm_dev = thick_membrane;				// membrane thickness, m
		double mffa_dev = mass_feedflow;			// mass flow rate @ feedair, kg/s
		double mfpf_dev = mass_permeateflow;		// mass flow rate @ permeateflow, kg/s
		int nl_dev = segment_length;				// segment number of L-direction
		int nd_dev = segment_deep;					// segment number of D-direction
		int nlay_dev = layernum_membrane;			// layer number of membrane
		int ftype_dev = flowtype_device;			// (number: 1,2,3)
		// ftype_dev: 1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or 3 or other: crossflow

		//MEMBRANE PROPERTY allocation
		double rp_membrane = poreradius_membrane;	// pore radius of membrane,m (the parameter will not be used if not selecting default model)
		double ps_membrane = porosity_membrane;		// porosity of membrane,(dimensionless) (the parameter will not be used if not selecting default model)
		double ts_membrane = tortuosity_membrane;	// tortuosity of membrane,(dimensionless) (the parameter will not be used if not selecting default model)
		double kt_membrane = thermocond_membrane;	// thermocondivity of membrane,kW/m-K (the parameter will not be used if not selecting default model)

		//FEEDFLOW AND PERMEATE FLOW allocation
		double t_fa = temp_feedflow;				// feedair temp, K
		double rh_fa = rh_feedflow;					// rh@feedair, (dimensionless)
		double hr_fa = airhumrat(t_fa, rh_fa);		// humid rate@fedair, kg-vap/kg-dry-air
		double vapfl_fa = mffa_dev / (1.0 + 1.0 / hr_fa);  // vapor flow rate @feedair, vap kg/s
		double p_fa = press_feedflow;				// feedair pressure, Pa
		double t_pf = temp_permeateflow;			// peameateflow temp, K
		double rh_pf = rh_permeateflow;				// rh@peameateflow, (dimensionless)
		double hr_pf = airhumrat(t_pf, rh_pf);		// humid rate@peameateflow, kg-vap/kg-dry-air
		double vapfl_pf = mfpf_dev / (1.0 + 1.0 / hr_pf);  // vapor flow rate @peameateflow, kg/s
		double p_pf = press_permeateflow;			// peameateflow pressure, Pa
		double enable_mtimpact = mt_impact_ht;		// 0:not enabled; 1: enabled mass permeated on ht impact

		int proppf_airvap = flowcondition_permeate;	// 0:vapor only in permeateflow; 1: air in permeateflow
		int propmem_model = membrane_model;			// 0:default model; 1: Dais data (rp_membrane,ps_membrane,ts_membrane,kt_membrane are invalue parameters); 2 constant value-model
		double htc_correct = htc_correction;		// HTC correction (0-1), (dimensionless)
		double ftc_correct = ftc_correction;		// FTC correction (0-1),(dimensionless)
		double memsurf_correct = memsurf_correction; //membrane surface correction for membrane surface deflection (0-inf, 1 means perfectly flat)
		if (proppf_airvap == 0) { t_pf = t_fa; }	// isothermal process considerd in the case of a vaccum permeate flow

		double * *map_temp_fa;
		double * *map_temp_fam;
		double * *map_flow_fa;
		double * *map_vap_fa;
		double * *map_humid_fa;
		double * *map_pvap_fa;
		double * *map_temp_pf;
		double * *map_temp_pfm;
		double * *map_flow_pf;
		double * *map_vap_pf;
		double * *map_humid_pf;
		double * *map_pvap_pf;
		double * *map_permeat_vap;
		if ((ftype_dev == 1) || (ftype_dev == 2))
		{
			map_temp_fa = dmatrix(0, 0, 0, nd_dev - 1);
			map_temp_fam = dmatrix(0, 0, 0, nd_dev - 1);
			map_flow_fa = dmatrix(0, 0, 0, nd_dev - 1);
			map_vap_fa = dmatrix(0, 0, 0, nd_dev - 1);
			map_humid_fa = dmatrix(0, 0, 0, nd_dev - 1);
			map_pvap_fa = dmatrix(0, 0, 0, nd_dev - 1);
			map_temp_pf = dmatrix(0, 0, 0, nd_dev - 1);
			map_temp_pfm = dmatrix(0, 0, 0, nd_dev - 1);
			map_flow_pf = dmatrix(0, 0, 0, nd_dev - 1);
			map_vap_pf = dmatrix(0, 0, 0, nd_dev - 1);
			map_humid_pf = dmatrix(0, 0, 0, nd_dev - 1);
			map_pvap_pf = dmatrix(0, 0, 0, nd_dev - 1);
			map_permeat_vap = dmatrix(0, 0, 0, nd_dev - 1);
		}
		else
		{
			map_temp_fa = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_temp_fam = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_flow_fa = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_vap_fa = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_humid_fa = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_pvap_fa = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_temp_pf = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_temp_pfm = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_flow_pf = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_vap_pf = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_humid_pf = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_pvap_pf = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
			map_permeat_vap = dmatrix(0, nl_dev - 1, 0, nd_dev - 1);
		}

		Component DaisHXM(l_dev, d_dev, hfa_dev, hpf_dev, tm_dev, mffa_dev, mfpf_dev, rp_membrane, ps_membrane, ts_membrane, kt_membrane,
			t_fa, vapfl_fa, p_fa, t_pf, vapfl_pf, p_pf, proppf_airvap, propmem_model, nl_dev, nd_dev, nlay_dev, ftype_dev, htc_correct, ftc_correct,
			map_temp_fa, map_temp_fam, map_vap_fa, map_humid_fa, map_flow_fa, map_pvap_fa, map_temp_pf, map_temp_pfm, map_vap_pf, map_humid_pf, map_flow_pf, map_pvap_pf, map_permeat_vap,
			enable_mtimpact, memsurf_correct);

		double tfa_in = t_fa;
		double tfa_out = DaisHXM.tfa_out;
		double humidfa_in = hr_fa;
		double humidfa_out = DaisHXM.vapfa_out / (mffa_dev - vapfl_fa);
		double rhfa_in = airrelhum(t_fa, hr_fa);
		double rhfa_out = airrelhum(DaisHXM.tfa_out, humidfa_out);
		double Pvapfa_in = vapPstaturate(t_fa)*airrelhum(t_fa, hr_fa);
		double Pvapfa_out = vapPstaturate(DaisHXM.tfa_out)*airrelhum(DaisHXM.tfa_out, humidfa_out);

		double tpf_in = t_pf;
		double tpf_out = DaisHXM.tpf_out;
		double humidpf_in = hr_pf;
		double humidpf_out = DaisHXM.vappf_out / (mfpf_dev - vapfl_pf);
		double rhpf_in = airrelhum(t_pf, hr_pf);
		double rhpf_out = airrelhum(DaisHXM.tpf_out, humidpf_out);
		double Pvappf_in = vapPstaturate(t_pf)*airrelhum(t_pf, hr_pf);
		double Pvappf_out = vapPstaturate(DaisHXM.tpf_out)*airrelhum(DaisHXM.tpf_out, humidpf_out);
		if (proppf_airvap == 0) 
		{humidpf_out = 1.0; rhpf_in = 1.0; rhpf_out = 1.0; Pvappf_in = p_pf; Pvappf_out = p_pf;}

		double mcpfa = aircp(tfa_in, humidfa_in)*mffa_dev;
		double mcppf = aircp(tpf_in, humidpf_in)*mfpf_dev;
		double effectiveness_l = mffa_dev*(humidpf_out-humidpf_in)/((humidfa_in-humidpf_in)*min(mffa_dev, mfpf_dev));
		double effectiveness_s = mcpfa*(tpf_out - tpf_in)/((tfa_in - tpf_in)*min(mcpfa, mcppf));
		if(proppf_airvap==0)
		{effectiveness_l = (Pvapfa_in - Pvapfa_out)/(Pvapfa_in- Pvappf_in); effectiveness_s = 0;}

		outputfile(nl_dev, nd_dev, nlay_dev, ftype_dev, tfa_in, tfa_out, humidfa_in, humidfa_out, rhfa_in, rhfa_out, Pvapfa_in, Pvapfa_out,
			tpf_in, tpf_out, humidpf_in, humidpf_out, rhpf_in, rhpf_out, Pvappf_in, Pvappf_out, map_temp_fa, map_temp_fam, map_vap_fa, map_humid_fa, map_flow_fa, map_pvap_fa,
			map_temp_pf, map_temp_pfm, map_vap_pf, map_humid_pf, map_flow_pf, map_pvap_pf, map_permeat_vap);

		cout << "tfa_in=" << tfa_in - 273.15 << "; humidfa_in=" << humidfa_in << "; rhfa_in=" << rhfa_in << "; Pvap_fain=" << Pvapfa_in << endl;
		cout << "tfa_out=" << tfa_out - 273.15 << " humidfa_out=" << humidfa_out << "; rhfa_out=" << rhfa_out << "; Pvap_faout=" << Pvapfa_out << endl;
		cout << "flfa_in=" << mffa_dev << "; flfa_out=" << DaisHXM.flfa_out << endl;
		cout << "tpf_in=" << tpf_in - 273.15 << "; humidpf_in=" << humidpf_in << "; rhpf_in=" << rhpf_in << "; Pvap_pfin=" << Pvappf_in << endl;
		cout << "tpf_out=" << tpf_out - 273.15 << " humidpf_out=" << humidpf_out << "; rhpf_out=" << rhpf_out << "; Pvap_pfout=" << Pvappf_out << endl;
		cout << "flpf_in=" << mfpf_dev << "; flpf_out=" << DaisHXM.flpf_out << endl;
		cout << "effectiveness_l=" << effectiveness_l << "; effectiveness_s=" << effectiveness_s << endl;

		if((ftype_dev==1)||(ftype_dev==2))
		{
			free_dmatrix(map_temp_fa, 0,0,0,nd_dev-1);
			free_dmatrix(map_temp_fam, 0,0,0,nd_dev-1);
			free_dmatrix(map_flow_fa, 0,0,0,nd_dev-1);
			free_dmatrix(map_vap_fa, 0,0,0,nd_dev-1);
			free_dmatrix(map_humid_fa, 0,0,0,nd_dev-1);
			free_dmatrix(map_pvap_fa, 0,0,0,nd_dev-1);
			free_dmatrix(map_temp_pf, 0,0,0,nd_dev-1);
			free_dmatrix(map_temp_pfm, 0,0,0,nd_dev-1);
			free_dmatrix(map_flow_pf, 0,0,0,nd_dev-1);
			free_dmatrix(map_vap_pf, 0,0,0,nd_dev-1);
			free_dmatrix(map_humid_pf, 0,0,0,nd_dev-1);
			free_dmatrix(map_pvap_pf, 0,0,0,nd_dev-1);
			free_dmatrix(map_permeat_vap, 0,0,0,nd_dev-1);
		}
		else
		{
			free_dmatrix(map_temp_fa, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_temp_fam, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_flow_fa, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_vap_fa, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_humid_fa, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_pvap_fa, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_temp_pf, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_temp_pfm, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_flow_pf, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_vap_pf, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_humid_pf, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_pvap_pf, 0,nl_dev-1,0,nd_dev-1);
			free_dmatrix(map_permeat_vap, 0,nl_dev-1,0,nd_dev-1);
		}
		system("pause");
	}
	catch(exception &e)
	{
		cout << "Caught Exception: " << e.what();
		system("pause");
		exit(1);
	}

	return 0;
}	// End of main
