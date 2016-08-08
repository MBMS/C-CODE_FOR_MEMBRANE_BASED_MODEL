//SegmentProperty.h written by Zhiming Gao (gaoz@ornl.gov)

#ifndef Segment_Property_Include
#define Segment_Property_Include

using namespace std;

struct SegmentProperty
{
	double m_rp, m_ps, m_ts, m_kt, m_temp;
	double fa_rho, fa_cp, fa_mu, fa_kt, fa_vap, fa_hr, fa_h, fa_temp, fa_flow, fa_htc, fa_mtc, fa_p, fa_vapdf;
	double pf_rho, pf_cp, pf_mu, pf_kt, pf_vap, pf_hr, pf_h, pf_temp, pf_flow, pf_htc, pf_mtc, pf_p, pf_vapdf;
	double wat_lhv, htc_correct, ftc_correct, memsurf_correct;
	int prop_airvap;

	SegmentProperty()
	{
		m_rp=0.25e-6;						// radius of membrane pore, (m)
		m_ps=0.55;							// membrane porosity, (dimensionless)
		m_ts=3.43;							// membrane tortuosity (dimensionless)
		m_kt=0.066e-3;						// membrane thermoconductivity (kW/m-K)
		m_temp=300.0;						// membrane temperature (K)
		fa_rho=0.9950;						// feed air density (kg/m3)
		fa_cp=1.007;						// feed air cp (kJ/kg-K)
		fa_mu=1.846e-4;						// feed air viscosity (Pa-s)
		fa_kt=0.0263e-3;					// feed air thermoconductivity (kW/m-K)
		fa_vap=0.0;							// feed air vap flow (kg/s)
		fa_hr=0.0;							// feed air humidity (kg/kg)
		fa_h=300.19;						// feed air enhalpy (kJ/kg)
		fa_temp=300.0;						// feed air temperature (K)
		fa_flow=1.0e-6;						// feed air mass flow per channel (Kg/s)
		fa_htc=1.0e-6;						// feed air heat transfer coefficient (kW/m2-K)
		fa_mtc=1.0e-6;						// feed air mass transfer coefficient (kg/m2-s)
		fa_p=1.01325e5;						// feed air pressure (Pa)
		fa_vapdf=2.82e-5;					// water vapor diffusion cofficent at feed side(m2/s)
		pf_rho=0.9950;						// permeate flow density (kg/m3)
		pf_cp=1.007;						// permeate flow  cp (kJ/kg-K)
		pf_mu=1.846e-4;						// permeate flow  viscosity (Pa-s)
		pf_kt=0.0263e-3;					// permeate flow  thermoconductivity (kW/m-K)
		pf_vap=0.0;							// permeate flow  vap flow (kg/s)
		pf_hr=0.0;							// permeate flow humidity (kg/kg)
		pf_h=300.19;						// permeate flow  enhalpy (kJ/kg)
		pf_temp=300;						// permeate flow  temperature (K)
		pf_flow=1.0e-6;						// permeate flow  mass flow per channel (Kg/s)
		pf_htc=1.0e-6;						// permeate flow  heat transfer coefficient (kW/m2-K)
		pf_mtc=1.0e-6;						// permeate flow  mass transfer coefficient (kg/m2-s)
		pf_p=1.01325e5;						// permeate flow  pressure (Pa)
		pf_vapdf=2.82e-5;					// water vapor diffusion cofficent at permeate side (m2/s)
		wat_lhv=2.438e3;					// water heat of vaporization (kJ/kg-K)

		htc_correct=1.0;					// HTC correction 
		ftc_correct=1.0;					// FTC correction 
		memsurf_correct = 1.0;				// membrane surface correction
		prop_airvap=1;						// property for permeate flow is moisture air (1) or vapor only (0)
	}

	SegmentProperty(const SegmentProperty &sgt)
	{
		m_rp=sgt.m_rp;
		m_ps=sgt.m_ps;
		m_ts=sgt.m_ts;
		m_kt=sgt.m_kt;
		fa_rho=sgt.fa_rho;	
		fa_cp=sgt.fa_cp;
		fa_mu=sgt.fa_mu;		
		fa_kt=sgt.fa_kt;	
		fa_vap=sgt.fa_vap;
		fa_hr=sgt.fa_hr;
		fa_h=sgt.fa_h;		
		fa_temp=sgt.fa_temp;
		fa_flow=sgt.fa_flow;
		fa_htc=sgt.fa_htc;	
		fa_mtc=sgt.fa_mtc;	
		fa_p=sgt.fa_p;	
		fa_vapdf=sgt.fa_vapdf;
		pf_rho=sgt.pf_rho;	
		pf_cp=sgt.pf_cp;		
		pf_mu=sgt.pf_mu;	
		pf_kt=sgt.pf_kt;	
		pf_vap=sgt.pf_vap;
		pf_hr=sgt.pf_hr;	
		pf_h=sgt.pf_h;		
		pf_temp=sgt.pf_temp;
		fa_flow=sgt.fa_flow;
		pf_htc=sgt.pf_htc;
		pf_mtc=sgt.pf_mtc;
		pf_p=sgt.pf_p;
		pf_vapdf=sgt.pf_vapdf;
		wat_lhv=sgt.wat_lhv;

		htc_correct = sgt.htc_correct;
		ftc_correct = sgt.ftc_correct;
		memsurf_correct = sgt.memsurf_correct;
		prop_airvap=sgt.prop_airvap;

	}

	SegmentProperty& operator=(const SegmentProperty &sgt)
	{
		if (this != &sgt)
		{
		m_rp=sgt.m_rp;
		m_ps=sgt.m_ps;
		m_ts=sgt.m_ts;
		m_kt=sgt.m_kt;
		fa_rho=sgt.fa_rho;	
		fa_cp=sgt.fa_cp;
		fa_mu=sgt.fa_mu;		
		fa_kt=sgt.fa_kt;	
		fa_vap=sgt.fa_vap;
		fa_hr=sgt.fa_hr;
		fa_h=sgt.fa_h;		
		fa_temp=sgt.fa_temp;
		fa_flow=sgt.fa_flow;
		fa_htc=sgt.fa_htc;	
		fa_mtc=sgt.fa_mtc;	
		fa_p=sgt.fa_p;
		fa_vapdf=sgt.fa_vapdf;
		pf_rho=sgt.pf_rho;	
		pf_cp=sgt.pf_cp;		
		pf_mu=sgt.pf_mu;	
		pf_kt=sgt.pf_kt;	
		pf_vap=sgt.pf_vap;	
		pf_hr=sgt.pf_hr;	
		pf_h=sgt.pf_h;		
		pf_temp=sgt.pf_temp;
		fa_flow=sgt.fa_flow;
		pf_htc=sgt.pf_htc;
		pf_mtc=sgt.pf_mtc;
		pf_vapdf=sgt.pf_vapdf;
		wat_lhv=sgt.wat_lhv;

		htc_correct = sgt.htc_correct;
		ftc_correct = sgt.ftc_correct;
		memsurf_correct = sgt.memsurf_correct;
		prop_airvap=sgt.prop_airvap;
		}
		return *this;
	}
};

namespace Constant_Coefficients
{
	const double pi_const=3.14159;		//dimensionless
	const double rgas_const=8.314;		//kJ/Kmol-K
	const double kilo_const=1.0e3;		//dimensionless
	const double mill_const=1.0e6;		//dimensionless
	const double molair=29.0;			//kg/kmol
	const double molwat=18.0;			//kg/kmol
}

#endif