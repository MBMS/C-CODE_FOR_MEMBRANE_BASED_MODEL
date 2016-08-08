// MoistAirVapProperty.cpp	Written by Zhiming Gao (gaoz@ornl.gov) 
// Modification and extension from a code obtained from Eckhard Groll at Purdue university, used with permission.
// Definition of namespace MoistAirVapProperty


#include <iostream>
#include <cmath>
#include <cstdlib>
#include "MoistAirVapProperty.h"

using namespace std;
using namespace Constant_Coefficients;


double MoistAirVapProperty::airh(double t, double rh)
{
	const double FS = 1.0039;
	const double P_ATM = 101325.0;
	
	if (rh < 0.0 || rh > 1.0)
	{
		cout << "Invalid relative humidity value (0-1): \t" << rh << endl;
		system("pause");
		exit(1);
	}

	if (t < 213.15 || t > 343.15)		// (t < 213.15 || t > 363.15)
	{
		cout << "Temperature out of range for function airh (213.15-343.15): " << t << endl;
		system("pause");
		exit(1);
	}
	
	double *c = new double[13];
	c[0] = -5.674359e3;
	c[1] = 6.3925247;
	c[2] = -9.6778430e-3;
	c[3] = 6.2215701e-7;
	c[4] = 2.0747825e-9;
	c[5] = -9.4840240e-13;
	c[6] = 4.1635019;
	c[7] = -5.8002206e3;
	c[8] = 1.3914993;
	c[9] = -4.8640239e-2;
	c[10] = 4.1764768e-5;
	c[11] = -1.4452093e-8;
	c[12] = 6.5459673;
	
	double pws = 0.0;
	if (t < 273.15)
	{pws = exp(c[0]/t + c[1] + c[2]*t + c[3]*pow(t, 2) + c[4]*pow(t, 3) + c[5]*pow(t, 4) + c[6]*log(t));}
	else
	{pws = exp(c[7]/t + c[8] + c[9]*t + c[10]*pow(t, 2) + c[11]*pow(t, 3) + c[12]*log(t));}
	
	double ws = 0.62198*((pws*FS)/(P_ATM - (FS*pws)));
	double degsat = rh/(1.0 + ((1.0 - rh)*ws)/0.62198);
	double w = degsat*ws;
	double enthalpy = (1.006*(t - 273.15) + w*(2501.0 + 1.805*(t - 273.15)));
	
	delete[] c;
	c = 0;
	
	return enthalpy;
}	// End of MoistAirVapProperty::airh


double MoistAirVapProperty::aircp(double t, double rh)
{
	const double FS = 1.0039;
	const double P_ATM = 101325.0;
	const double DELTA_T = 0.1;
	
	if (rh < 0.0 || rh > 1.0)
	{
		cout << "Invalid relative humidity value(0-1): \t" << rh << endl;
		system("pause");
		exit(1);
	}

	if (t < 213.15 || t > 343.15)
	{
		cout << "Temperature out of range for function aircp (213.15-343.15): " << t << endl;
		system("pause");
		exit(1);
	}
	
	double *c = new double[13];
	c[0] = -5.674359e3;
	c[1] = 6.3925247;
	c[2] = -9.6778430e-3;
	c[3] = 6.2215701e-7;
	c[4] = 2.0747825e-9;
	c[5] = -9.4840240e-13;
	c[6] = 4.1635019;
	c[7] = -5.8002206e3;
	c[8] = 1.3914993;
	c[9] = -4.8640239e-2;
	c[10] = 4.1764768e-5;
	c[11] = -1.4452093e-8;
	c[12] = 6.5459673;
	
	double pws = 0.0;
	if (t < 273.15)
	{pws = exp(c[0]/t + c[1] + c[2]*t + c[3]*pow(t, 2) + c[4]*pow(t, 3) + c[5]*pow(t, 4) + c[6]*log(t));}
	else
	{pws = exp(c[7]/t + c[8] + c[9]*t + c[10]*pow(t, 2) + c[11]*pow(t, 3) + c[12]*log(t));}
	
	double ws = 0.62198*((pws*FS)/(P_ATM - (FS*pws)));
	double degsat = rh/(1.0 + ((1.0 - rh)*ws)/0.62198);
	double w = degsat*ws;
	double tlow = t - DELTA_T;
	double thigh = t + DELTA_T;
	double hlow = 1.006*(tlow - 273.15) + w*(2501.0 + 1.805*(tlow - 273.15));
	double hhigh = 1.006*(thigh - 273.15) + w*(2501.0 + 1.805*(thigh - 273.15));
	double specheat = ((hhigh - hlow)/(thigh - tlow));
	
	delete[] c;
	c = 0;
	
	return specheat;
}	// End of MoistAirVapProperty::aircp


double MoistAirVapProperty::airrelhum(double t, double w)
{
	const double FS = 1.0039;
	const double P_ATM = 101325.0;
	
	double *c = new double[13];
	c[0] = -5.674359e3;
	c[1] = 6.3925247;
	c[2] = -9.6778430e-3;
	c[3] = 6.2215701e-7;
	c[4] = 2.0747825e-9;
	c[5] = -9.4840240e-13;
	c[6] = 4.1635019;
	c[7] = -5.8002206e3;
	c[8] = 1.3914993;
	c[9] = -4.8640239e-2;
	c[10] = 4.1764768e-5;
	c[11] = -1.4452093e-8;
	c[12] = 6.5459673;
	
	double pws = 0.0;
	if (t < 273.15)
	{pws = exp(c[0]/t + c[1] + c[2]*t + c[3]*pow(t, 2) + c[4]*pow(t, 3) + c[5]*pow(t, 4) + c[6]*log(t));}
	else
	{pws = exp(c[7]/t + c[8] + c[9]*t + c[10]*pow(t, 2) + c[11]*pow(t, 3) + c[12]*log(t));}
	
	double ws = 0.62198*((pws*FS)/(P_ATM - (FS*pws)));
	double degsat = w/ws;
	double rh = degsat/(1.0 - (1.0 - degsat)*(FS*pws/P_ATM));
	if (rh > 1.0)
	{rh = 1.0;}
	
	delete[] c;
	c = 0;
	
	return rh;	
}	// End of MoistAirVapProperty::relhum


double MoistAirVapProperty::airhumrat(double t, double rh)
{
	const double FS = 1.0039;
	const double P_ATM = 101325.0;
	
	if (rh < 0.0 || rh > 1.0)
	{
		cout << "Invalid relative humidity value(0-1): \t" << rh << endl;
		system("pause");
		exit(1);
	}

	if (t < 213.15 || t > 343.15)
	{
		cout << "Temperature out of range (213.15-343.15) for function humrat: " << t << endl;
		system("pause");
		exit(1);
	}
	
	double *c = new double[13];
	c[0] = -5.674359e3;
	c[1] = 6.3925247;
	c[2] = -9.6778430e-3;
	c[3] = 6.2215701e-7;
	c[4] = 2.0747825e-9;
	c[5] = -9.4840240e-13;
	c[6] = 4.1635019;
	c[7] = -5.8002206e3;
	c[8] = 1.3914993;
	c[9] = -4.8640239e-2;
	c[10] = 4.1764768e-5;
	c[11] = -1.4452093e-8;
	c[12] = 6.5459673;
	
	double pws = 0.0;
	if (t < 273.15)
	{pws = exp(c[0]/t + c[1] + c[2]*t + c[3]*pow(t, 2) + c[4]*pow(t, 3) + c[5]*pow(t, 4) + c[6]*log(t));}
	else
	{pws = exp(c[7]/t + c[8] + c[9]*t + c[10]*pow(t, 2) + c[11]*pow(t, 3) + c[12]*log(t));}
	
	double ws = 0.62198*((pws*FS)/(P_ATM - (FS*pws)));
	double degsat = rh/(1.0 + ((1.0 - rh)*ws)/0.62198);
	double w = degsat*ws;
	
	delete[] c;
	c = 0;
	
	return w;
}	// End of MoistAirVapProperty::humrat


double MoistAirVapProperty::airmu(double t)
{
	double mu = 4.49060743178e-14*pow(t, 3) - 7.7417957043e-11*pow(t, 2) + 8.131519014305697e-8*t - 2.095906817700090e-7;
	return mu;
}	// End of MoistAirVapProperty::airmu


double MoistAirVapProperty::airk(double t)
{
	double k = 1.348896201837e-10*pow(t, 3) - 1.5925643107e-7*pow(t, 2) + 1.2892406415e-4*t	- 1.843581337921997e-3;
	k *= 1.0e-3;	// Convert k from W/mk to kW/mk
	return k;
}	// End of MoistAirVapProperty::airk


double MoistAirVapProperty::airpr(double t, double rh)
{
	double pr = aircp(t, rh)*airmu(t)/airk(t);
	return pr;
}	// End of MoistAirVapProperty::airpr


double MoistAirVapProperty::airthsat(double h_sat)
{
	double t = 0.0;
	if (h_sat < -60.334 || h_sat > 3867.599)	// (h_sat < -35.0 || h_sat > 170.0)
	{
		cout << "Enthalpy out of range (-60.334~3867.599) for function airthsat: " << h_sat << endl;
		system("pause");
		exit(1);
	}
	else
	{
		t = -1.4797e-11*pow(h_sat, 4) + 1.1899e-7*pow(h_sat, 3) - 3.1245e-4*pow(h_sat, 2) + 0.3062*h_sat + 259.4960; // Ortiz fit
	}
	return t;
}	// End of MoistAirVapProperty::airthsat


double MoistAirVapProperty::airthdry(double h_dry)
{
	double t = 0.0;
	if (h_dry < -60.351 || h_dry > 90.681)
	{
		cout << "Enthalpy out of range (-60.351~90.681) for function airthdry." << h_dry << endl;
		system("pause");
		exit(1);
	}
	else
	{
		t = -3.1907e-11*pow(h_dry, 4) - 1.2588e-7*pow(h_dry, 3) - 7.2107e-6*pow(h_dry, 2) + 0.9942*h_dry + 273.1499;
	}
	return t;
}	// End of MoistAirVapProperty::airthdry


double MoistAirVapProperty::airrho(double t)
{
	double rho = -2.000119978063e-9*pow(t, 4) + 2.337924328125e-6*pow(t, 3) - 1.0039106594e-3*pow(t, 2) + 1.831951129154464e-1*t - 1.035534408474814e1;
	return rho;
}	// End of MoistAirVapProperty::airrho


double MoistAirVapProperty::airthrh(double h, double rh)
{
	const double T_TOL = 1.0e-1;
	const double BIG_NEGATIVE_NUMBER = -1.11e30;
		
	double t_1 = airthsat(h);
	double t_2 = airthdry(h);
	double t = 0.0;
	double f = 0.0;
	double lastGuess = BIG_NEGATIVE_NUMBER;
	
	while (fabs(t - lastGuess) > T_TOL)
	{
	 	lastGuess = t;
		t = (t_1 + t_2)/2.0;
		f = h - airh(t, rh);
	 	if (f < 0.0)
	 	{t_2 = t;}
	 	else
	 	{t_1 = t;}
	}
	return t;
}	// End of MoistAirVapProperty::airthrh


double MoistAirVapProperty::airthw(double h, double w)
{
	const double T_TOL = 1.0e-1;
	const double BIG_NEGATIVE_NUMBER = -1.11e30;
		
	double t_1 = airthsat(h);
	double t_2 = airthdry(h);
	double t = 0.0;
	double f = 0.0;
	double lastGuess = BIG_NEGATIVE_NUMBER;
	
	while (fabs(t - lastGuess) > T_TOL)
	{
	 	lastGuess = t;
		t = (t_1 + t_2)/2.0;
		double rh = airrelhum(t, w);
		f = h - airh(t, rh);
	 	if (f < 0.0)
	 	{t_2 = t;}
	 	else
	 	{t_1 = t;}
	}
	return t;
}	// End of MoistAirVapProperty::airthw


double MoistAirVapProperty::airrhth(double t, double h)
{
	const double RH_TOL = 1.0e-2;
	const double BIG_NEGATIVE_NUMBER = -1.11e30;
		
	double rh_1 = 0.0;
	double rh_2 = 1.0;
	double rh = 0.0;
	double f = 0.0;
	double lastGuess = BIG_NEGATIVE_NUMBER;
	
	while (fabs(rh - lastGuess) > RH_TOL)
	{
	 	lastGuess = rh;
		rh = (rh_1 + rh_2)/2.0;
		f = h - airh(t, rh);
	 	if (f < 0.0)
	 	{rh_2 = rh;}
	 	else
	 	{rh_1 = rh;}
	}
	return rh;
}	// End of MoistAirVapProperty::airrhth


double MoistAirVapProperty::airgetcs(double t_sat)
{
	const double DELTA_T_SAT = 0.1;
	
	double rhsat = 1.0;
	double cstlow = t_sat - DELTA_T_SAT;
	double csthi = t_sat + DELTA_T_SAT;
	double cshlow = airh(cstlow, rhsat);
	double cshhi = airh(csthi, rhsat);
	double cpsat = ((cshhi - cshlow)/(csthi - cstlow));
	
	return cpsat;
}	// End of MoistAirVapProperty::getcs


double MoistAirVapProperty::airtdp(double t, double rh)
{
	const double T_TOL = 1.0e-1;
	const double DELTA_T = 10.0;
	const double HIGH_AIR_RH_LIMIT = 1.0;
	const double BIG_NEGATIVE_NUMBER = -1.11e30;
	
	double w = airhumrat(t, rh);
	double t_1 = t;
	double t_2 = t - DELTA_T;
	double f_1 = w - airhumrat(t_1, HIGH_AIR_RH_LIMIT);
	double f_2 = w - airhumrat(t_2, HIGH_AIR_RH_LIMIT);
	double dewpt = 0.0;
	double f = 0.0;
	double lastGuess = BIG_NEGATIVE_NUMBER;
	
	while (fabs(dewpt - lastGuess) > T_TOL)
	{
	 	lastGuess = dewpt;
		dewpt = (t_1 + t_2)/2.0;
	 	f = f_1 - ((f_2 - f_1)/(t_2 - t_1))*(t_1 - dewpt);
	 	if (f > 0.0)
	 	{t_1 = dewpt;}
	 	else
	 	{t_2 = dewpt;}
	}
	return dewpt;
}	// End of MoistAirVapProperty::tdp


double MoistAirVapProperty::airdiff(double t)
{
	double diff_vapair=2.82e-3;
	return diff_vapair;
}	//End of MoistAirVapProperty::airdiff


double MoistAirVapProperty::airre_flatplate(double t, HXM* hExchanger, int flowside)
{
	//flowside: 1-> feed airflow side; 2-> permeate flow side;
	if(flowside == 1)
	{
		double mdot=hExchanger->get_mflwseg_fa();
		double a_face = hExchanger->aseg_crossfa();
		double d_c =hExchanger->get_lcseg_fa();
		double mu = airmu(t);
	
		double gpm2_air = mdot/(a_face);
		double re = gpm2_air*d_c /mu;
		return re;
	}
	else if(flowside == 2)
	{
		double mdot=hExchanger->get_mflwseg_pf();
		double a_face = hExchanger->aseg_crosspf();
		double d_c =hExchanger->get_lcseg_pf();
		double mu = airmu(t);
	
		double gpm2_air = mdot/(a_face);
		double re = gpm2_air*d_c /mu;
		return re;
	}
	else
	{
		cout << "Invalid flowside value (1: feed airflow; 2: permeate flow:): \t" << flowside << endl;
		system("pause");
		exit(1);
	}
}	// End of MoistAirVapProperty::airre_flatplate


double MoistAirVapProperty::airffr_flatplate(SegmentProperty &segment, HXM *hExchanger, int flowside)
{
	double f = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);
	if (flowside == 1)
	{
		//feed side;
		double temp=segpoint.fa_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		if(re<2300)
		{f=96.0/re;}
		else if (re<2.0e4)
		{f=0.316*pow(re, -0.25);}
		else
		{f=0.184*pow(re, -0.20);}

		return f*segpoint.ftc_correct;
	}
	else if (flowside == 2)
	{
		//permeate side;
		double temp=segpoint.fa_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		if(re<2300)
		{f=96.0/re;}
		else if (re<2.0e4)
		{f=0.316*pow(re, -0.25);}
		else
		{f=0.184*pow(re, -0.20);}

		return f*segpoint.ftc_correct;
	}
	else
	{
		cout << "Invalid flowside value (1: feed airflow; 2: permeate flow:): \t" << flowside << endl;
		system("pause");
		exit(1);
	}
}	// End of MoistAirVapProperty::ffr_flatplate


double MoistAirVapProperty::airhtc_flatplate(SegmentProperty &segment, HXM *hExchanger, int flowside)
{
	double htc = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);
	if (flowside == 1)
	{
		//feed side;
		double temp=segpoint.fa_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		double rh = airrelhum(temp, segpoint.fa_hr);
		double pr = airpr(temp, rh);
		double k  =airk(temp);
		double lc =hExchanger->get_lcseg_fa();

		if(re<2300)
		{htc= 7.54*k/lc;}
		else
		{htc = 0.023*k*pow(re, 0.8)*pow(pr, 0.3)/lc;}
		return htc*segpoint.htc_correct;		
	}
	else if (flowside == 2)
	{
		//permeate side;
		double temp=segpoint.pf_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		double rh = airrelhum(temp, segpoint.pf_hr);
		double pr = airpr(temp, rh);
		double k  = airk(temp);
		double lc =hExchanger->get_lcseg_pf();
		
		if(re<2300)
		{htc= 7.54*k/lc;}
		else
		{htc = 0.023*k*pow(re, 0.8)*pow(pr, 0.4)/lc;}
		return htc*segpoint.htc_correct;		
	}
	else
	{
		cout << "Invalid flowside value (1: feed airflow; 2: permeate flow:): \t" << flowside << endl;
		system("pause");
		exit(1);
	}
}	// End of MoistAirVapProperty::airhtc_flatplate


double MoistAirVapProperty::airmtc_flatplate(SegmentProperty &segment, HXM *hExchanger, int flowside)
{
	double mtc = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);
	if (flowside == 1)
	{
		//feed side;
		double temp=segpoint.fa_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		//double rh = relhum(temp, segpoint.fa_hr);
		double df = airdiff(temp);
		double sc = airmu(temp)/df;
		double lc = hExchanger->get_lcseg_fa();

		if(re<2300)
		{mtc=7.54*df/lc;}
		else
		{mtc = 0.023*df*pow(re, 0.8)*pow(sc, 0.4)/lc;}
		return mtc*segpoint.htc_correct;
	}
	else if (flowside == 2)
	{
		//permeate side;
		double temp=segpoint.pf_temp;
		double re = airre_flatplate(temp, hExchanger, flowside);
		//double rh = relhum(temp, segpoint.pf_hr);
		double df = airdiff(temp);
		double sc = airmu(temp)/df;
		double lc =hExchanger->get_lcseg_pf();
		
		if(re<2300)
		{mtc=7.54*df/lc;}
		else
		{mtc = 0.023*df*pow(re, 0.8)*pow(sc, 0.4)/lc;}
		return mtc*segpoint.htc_correct;		
	}
	else
	{
		cout << "Invalid flowside value (1: feed airflow; 2: permeate flow:): \t" << flowside << endl;
		system("pause");
		exit(1);
	}
}	// End of MoistAirVapProperty::airmtc_flatplate


void MoistAirVapProperty::allprop_flatplate(SegmentProperty &segment, HXM *hExchanger, int flowside)
{
	SegmentProperty segpoint = SegmentProperty(segment);
	if (flowside == 1)
	{
		//feed side;
		double temp=segpoint.fa_temp;
		double rh = airrelhum(temp, segpoint.fa_hr);
		segment.fa_rho=airrho(temp);
		segment.fa_cp=aircp(temp, rh);
		segment.fa_mu= airmu(temp);	
		segment.fa_kt=airk(temp);
		segment.fa_vapdf=airdiff(temp);
		segment.fa_h=airh(temp, rh);
		segment.fa_htc= airhtc_flatplate(segment, hExchanger, flowside);
		segment.fa_mtc= airmtc_flatplate(segment, hExchanger, flowside);
	}
	else if (flowside == 2)
	{
		//permeate side;
		double temp = segpoint.pf_temp;
		if (segment.prop_airvap == 0) 
		{
			segment.pf_rho = vaprho(temp);
			segment.pf_cp = vapcp(temp);
			segment.pf_mu = vapmu(temp);
			segment.pf_kt = vapk(temp);
			segment.pf_vapdf = vapdiff(temp);
			segment.pf_h = vaph(temp);
			segment.pf_htc = vaphtc_flatplate(segment, hExchanger);
			segment.pf_mtc = vapmtc_flatplate(segment, hExchanger);
		}
		else
		{
			double rh = airrelhum(temp, segpoint.pf_hr);
			segment.pf_rho = airrho(temp);
			segment.pf_cp = aircp(temp, rh);
			segment.pf_mu = airmu(temp);
			segment.pf_kt = airk(temp);
			segment.pf_vapdf = airdiff(temp);
			segment.pf_h = airh(temp, rh);
			segment.pf_htc = airhtc_flatplate(segment, hExchanger, flowside);
			segment.pf_mtc = airmtc_flatplate(segment, hExchanger, flowside);
		}
	}
	else
	{
		cout << "Invalid flowside value (1: feed airflow; 2: permeate flow:): \t" << flowside << endl;
		system("pause");
		exit(1);
	}
} // End of MoistAirVapProperty::allprop_flatplate


double MoistAirVapProperty::vaph(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	double h_vap = vapcp(t) * t;
	return h_vap;
} // End of MoistAirVapProperty::vaph


double MoistAirVapProperty::vapcp(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	//double cp_vap = 1.864;
	double cp_vap = 8.8889E-8*pow(t, 3) - 7.3905E-5*pow(t, 2) + 2.1568E-2*t - 3.0484E-1;
	return cp_vap;
} // End of MoistAirVapProperty::vapcp


double MoistAirVapProperty::vapmu(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	double mu_vap = -2.6759E-13*pow(t, 3) + 2.8672E-10*pow(t, 2) - 6.8010E-8*t + 1.1582E-5;
	return mu_vap;
} // End of MoistAirVapProperty::vapmu


double MoistAirVapProperty::vapk(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	double k_vap = 3.5185E-13*pow(t, 3) - 1.9266E-10*pow(t, 2) + 9.0888E-8*t - 8.6417E-7;
	return k_vap;
} // End of MoistAirVapProperty::vapk


double MoistAirVapProperty::vappr(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	double pr = vapcp(t) * vapmu(t) / vapk(t);
	return pr;
}	// End of MoistAirVapProperty::vappr


double MoistAirVapProperty::vaprho(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}	
	double rho_vap = -11.201 + 1.2016e-1 * t - 4.3151e-4 * pow(t, 2) + 5.1906e-7 * pow(t, 3);
	return rho_vap;
} 	// End of MoistAirVapProperty::vapho


double MoistAirVapProperty::vapdiff(double t)
{
	if(t<280.0){t=280.0;}
	if(t>330.0){t=330.0;}
	double diff_vap = 2.180733e4 - 2.719014e2 * t + 1.273726 * pow(t, 2) - 2.656102e-3 * pow(t, 3) + 2.079792e-6 * pow(t, 4);
	return diff_vap/10000.0;
} 	// End of MoistAirVapProperty::vapdiff


double MoistAirVapProperty::vapPstaturate(double t)
{
	double ps_f = exp(77.3450 + 0.0057 * t - 7236 / t)/pow(t,8.2);
	return ps_f;
}	// End of MoistAirVapProperty::waterViscosity


double MoistAirVapProperty::waterLatentHeat(double t)
{
	double h_fg = 2507.5489 - 2.553829*(t - 273.15);
	return h_fg;
}	// End of MoistAirVapProperty::waterLatentHeat


double MoistAirVapProperty::vapre_flatplate(double t, HXM* hExchanger)
{
	double mdot = hExchanger->get_mflwseg_pf();
	double a_face = hExchanger->aseg_crosspf();
	double d_c = hExchanger->get_lcseg_pf();
	double mu = vapmu(t);

	double gpm2_air = mdot / (a_face);
	double re = gpm2_air*d_c / mu;
	return re;	
}	// End of MoistAirVapProperty::vapre_flatplate


double MoistAirVapProperty::vapffr_flatplate(SegmentProperty &segment, HXM *hExchanger)
{
	double f = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);

	//permeate side;
	double temp = segpoint.fa_temp;
	double re = vapre_flatplate(temp, hExchanger);
	if (re<2300)
	{f = 96.0 / re;}
	else if (re<2.0e4)
	{f = 0.316*pow(re, -0.25);}
	else
	{f = 0.184*pow(re, -0.20);}
	return f;
}	// End of MoistAirVapProperty::ffr_flatplate


double MoistAirVapProperty::vaphtc_flatplate(SegmentProperty &segment, HXM *hExchanger)
{
	double htc = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);
	
	//permeate side;
	double temp = segpoint.pf_temp;
	double re = vapre_flatplate(temp, hExchanger);
	double pr = vappr(temp);
	double k = vapk(temp);
	double lc = hExchanger->get_lcseg_pf();

	if (re<2300)
	{htc = 7.54*k / lc;}
	else
	{htc = 0.023*k*pow(re, 0.8)*pow(pr, 0.4) / lc;}
	return htc;
}	// End of MoistAirVapProperty::htc_flatplate


double MoistAirVapProperty::vapmtc_flatplate(SegmentProperty &segment, HXM *hExchanger)
{
	double mtc = 0.0;
	SegmentProperty segpoint = SegmentProperty(segment);
	
	//permeate side;
	double temp = segpoint.pf_temp;
	double re = vapre_flatplate(temp, hExchanger);
	double df = vapdiff(temp);
	double sc = vapmu(temp) / df;
	double lc = hExchanger->get_lcseg_pf();
	
	if (re<2300)
	{mtc = 7.54*df / lc;}
	else
	{mtc = 0.023*df*pow(re, 0.8)*pow(sc, 0.4) / lc;}
	return mtc;
}	// End of MoistAirVapProperty::mtc_flatplate
