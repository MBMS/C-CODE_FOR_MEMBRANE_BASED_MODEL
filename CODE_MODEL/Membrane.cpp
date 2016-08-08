// Membrane.cpp	Written by Zhiming Gao (gaoz@ornl.gov).
// Definition of namespace Membrane


#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Membrane.h"


using namespace std;
using namespace Constant_Coefficients;
using namespace MoistAirVapProperty;


double Membrane::knudsen_diffus(SegmentProperty &segment, HXM* hExchanger, double t_membrane)
{
	SegmentProperty segpoint = SegmentProperty(segment);
	double dmembrane = hExchanger->thick_membrane;
	double c_kd1=segpoint.m_rp*segpoint.m_ps/(segpoint.m_ts*dmembrane);
	double c_kd2=molwat/(rgas_const*t_membrane)/1000.0;
	double c_kd =1.064*c_kd1*pow(c_kd2, 0.5);

	return c_kd;
}	// End of Membrane::knudsen_diffus

double Membrane::molecular_diffus(SegmentProperty &segment, HXM* hExchanger, double t_membrane, double p_mvap, double p_mair)
{
	SegmentProperty segpoint = SegmentProperty(segment);
	double dmembrane = hExchanger->thick_membrane;
	double c_md1=segpoint.m_ps/(segpoint.m_ts*dmembrane);
	double c_md2=21.2e-6*(1.0+0.0071*(t_membrane-273.15))*p_mvap/p_mair;
	double c_md3=molwat/(rgas_const*t_membrane)/1000.0;
	double c_md =c_md1*c_md2*c_md3;

	return c_md;
}	// End of Membrane::knudsen_diffus

double Membrane::poise_flow(SegmentProperty &segment, HXM* hExchanger, double t_membrane, double p_mvap)
{
	SegmentProperty segpoint = SegmentProperty(segment);
	double dmembrane = hExchanger->thick_membrane;
	double mu_mvap = vapmu(t_membrane);
	double c_pf1=pow(segpoint.m_rp, 2.0)*segpoint.m_ps/(segpoint.m_ts*dmembrane);
	double c_pf2=molwat*p_mvap/(mu_mvap*rgas_const*t_membrane)/1000.0;
	double c_pf =0.125*c_pf1*c_pf2;

	return c_pf;
}	// End of Membrane::knudsen_diffus

double Membrane::watpermeatd(SegmentProperty &segment, HXM* hExchanger, double p_fa, double p_pf, double t_membrane, double p_mvap, double p_mair)
{
	double c_kd = knudsen_diffus(segment, hExchanger,  t_membrane);
	double c_md = molecular_diffus(segment, hExchanger,  t_membrane, p_mvap, p_mair);
	double c_pf = poise_flow(segment, hExchanger,  t_membrane,  p_mvap);
	double ck = segment.m_rp*segment.m_ps / (segment.m_ts*hExchanger->thick_membrane);
	double cm = segment.m_ps / (segment.m_ts*hExchanger->thick_membrane);
	double cp = segment.m_rp*segment.m_rp*segment.m_ps / (segment.m_ts*hExchanger->thick_membrane);
	double kvwat=(c_kd*c_md/(c_kd+c_md)+c_pf);	//kg/m2-s-pa
	double fvwat=kvwat*(p_fa - p_pf);			//Based on vapor pressure drive force
	if (p_fa - p_pf < 0) fvwat = 0;
	return fvwat;
}	// End of Membrane::watpermeat

double Membrane::watpermeatc(double p_fa, double p_pf)
{
	double fvwat= 8.0e-7*(p_fa - p_pf);		//Based on vapor pressure drive force
	if (p_fa - p_pf < 0) fvwat = 0;
	return fvwat;
}	// End of Membrane::watpermeatc @ given constant permeation specific flux

double Membrane::watpermeatdais(double p_fa, double p_pf, double rh_fa, double fl_fa, double t_fa, double d_fa)
{

	double t25=273.15+25;
	double t30=273.15+30;
	double t35=273.15+35;
	double j25_40 = 2.48e-4*pow((rh_fa*100.0), 1.71);
	double j30_40 = 2.48e-4*pow((rh_fa*100.0), 1.73);
	double j35_40 = 1.76e-4*pow((rh_fa*100.0), 1.73);
	double j25_100 = 2.48e-4*pow((rh_fa*100.0), 1.76);
	double j30_100 = 2.53e-4*pow((rh_fa*100.0), 1.78);
	double j35_100 = 1.24e-4*pow((rh_fa*100.0), 1.90);
	double jt = 0.0;

	double re_fa = fl_fa*d_fa/airmu(t_fa);
	double re_fa30_40 = airrho(t30)*40/6.0e4/4.57e-4*2.0e-3/airmu(t30);
	double re_fa30_100= airrho(t30)*100/6.0e4/4.57e-4*2.0e-3/airmu(t30);

	if(t_fa<t30)
	{
		double dt=(t_fa-t25)/(t30-t25);
		double jt40=(j30_40-j25_40)*dt+j25_40;
		double jt100=(j30_100-j25_100)*dt+j25_100;
		//if (re_fa > re_fa30_100) { re_fa = re_fa30_100; }
		//if (re_fa < re_fa30_40) { re_fa = re_fa30_40; }
		jt=jt40+(jt100-jt40)*(re_fa-re_fa30_40)/(re_fa30_100-re_fa30_40);	//kg/hr-m2-kpa
	}
	else 
	{
		double dt=(t_fa-t30)/(t35-t30);
		double jt40=(j35_40-j30_40)*dt+j30_40;
		double jt100=(j35_100-j30_100)*dt+j30_100;
		//if (re_fa > re_fa30_100) { re_fa = re_fa30_100; }
		//if (re_fa < re_fa30_40) { re_fa = re_fa30_40; }
		jt=jt40+(jt100-jt40)*(re_fa-re_fa30_40)/(re_fa30_100-re_fa30_40);	//kg/hr-m2-kpa
	}
	double fvwat = jt*(p_fa - p_pf)/3.6e6;		//Based on vapor pressure drive force
	if (p_fa - p_pf < 0) fvwat = 0;
	return fvwat;
}	// End of Membrane::watpermeatdais

double Membrane::airpermeat(double p_fa, double p_pf, double w_fa, double w_pf)
{
	double kvair=0.0;		//kg/m2-s-pa
	double fvair1=kvair;	//shall remove molair???
	double fvair2=p_fa*molwat*(w_fa/(1.0+w_fa))/(molwat/(1.0+w_fa)+molair*w_fa/(1.0+w_fa));
	double fvair3=p_pf*molwat*(w_pf/(1.0+w_pf))/(molwat/(1.0+w_pf)+molair*w_pf/(1.0+w_pf));
	double fvair=fvair1*(fvair2-fvair3);
	return fvair;
}	// End of Membrane::airpermeat