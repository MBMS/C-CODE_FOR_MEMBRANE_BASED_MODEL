//Component.cpp written by Zhiming Gao (gaoz@ornl.gov)


#include <iostream>
#include <cstdlib>
#include "Component.h"

using namespace std;
using namespace Constant_Coefficients;
using namespace MoistAirVapProperty;
using namespace Membrane;


Component::Component(double lvalue, double dvalue, double hfavalue, double hpfvalue, double tmvalue, double mffa, double mfpf, 
		  double rpmvalue, double psmvalue, double tsmvalue, double ktmvalue, 
		  double tfa_in, double vapfa_in, double pfa_in, double tpf_in, double vappf_in, double ppf_in,
		  int prop_pfav, int prop_memmodel, int lnvalue, int dnvalue, int laynvalue, int ftype, double htc_correct, double ftc_correct,
		  double * *map_temp_fa, double * *map_temp_fam, double * *map_vap_fa, double * *map_humid_fa, double * *map_flow_fa, double * *map_pvap_fa,
		  double * *map_temp_pf, double * *map_temp_pfm, double * *map_vap_pf, double * *map_humid_pf, double * *map_flow_pf, double * *map_pvap_pf, double * *map_permeat_vap,
		  double enable_mtimpact, double memsurf_correct)
{
	length=lvalue;
	depth =dvalue;
	thick_fa=hfavalue; 
	thick_pf=hpfvalue;
	thick_membrane=tmvalue;
	mf_fa=mffa;
	mfvap_fa=vapfa_in;
	mf_pf=mfpf;
	mfvap_pf=vappf_in;

	rp_membrane=rpmvalue;
	ps_membrane=psmvalue;
	ts_membrane=tsmvalue;	
	kt_membrane=ktmvalue;

	num_lsegment=lnvalue;
	num_dsegment=dnvalue;
	num_layer=laynvalue;
	proppf_airvap=prop_pfav;
	htc_correction= htc_correct;
	ftc_correction= ftc_correct;
	memsurf_correction = memsurf_correct;

	flowtype=ftype;			// 1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or other: crossflow
	int numchannel_fa = num_layer/2+1;
	int numchannel_pf = (num_layer+1)/2;
	if (num_layer>2) {numchannel_fa=(num_layer+1)/2; numchannel_pf=(num_layer+1)/2;}

	HXM HXComponent(length, depth, thick_fa, thick_pf, thick_membrane, mf_fa, mfvap_fa, mf_pf, mfvap_pf, memsurf_correction, num_lsegment, num_dsegment, num_layer, flowtype);
	HXM* Dehumidifier=  &HXComponent;
	double mflwseg_fa_value=Dehumidifier->get_mflwseg_fa();
	double mflwvapseg_fa_value=Dehumidifier->get_mflwvapseg_fa();
	double lseg_fa_value=Dehumidifier->get_lseg_fa();
	double dseg_fa_value=Dehumidifier->get_dseg_fa();
	double mflwseg_pf_value=Dehumidifier->get_mflwseg_pf();	
	double mflwvapseg_pf_value=Dehumidifier->get_mflwvapseg_pf();	
	double lseg_pf_value=Dehumidifier->get_lseg_pf();
	double dseg_pf_value=Dehumidifier->get_dseg_pf();
	double delta_mem_value=thick_membrane;

	SegmentProperty segpoint;
	if (flowtype == 1)
	{
		//             fa_in|pf_in      : parallel uncounterflow
		//			 ______|____|____
		//		Hpf	|      |    v   |\ 
		//		Hm	|======v========| \
		//		Hfa	|_______________|\ \
		//			\      L        \D\|
		//			 \           |   \ |
		//			  \__________|____\|
		//				   |     V pf_out
		//				   |
		//				   V fa_out
		// 
		//cout<<"Notice: it is a parallel uncounterflow along D dirction"<<endl; 
		temp_fa=tfa_in;
		flow_fa= mflwseg_fa_value;
		vap_fa= mflwvapseg_fa_value;
		humid_fa=mflwvapseg_fa_value/(mflwseg_fa_value-mflwvapseg_fa_value);
		press_fa=pfa_in; 
		temp_pf=tpf_in;
		flow_pf = mflwseg_pf_value;
		vap_pf=mflwvapseg_pf_value;
		if(proppf_airvap == 0) {humid_pf=1.0;}
		if(proppf_airvap == 1) {humid_pf=mflwvapseg_pf_value/(mflwseg_pf_value-mflwvapseg_pf_value);}
		press_pf=ppf_in;

		for(int i=0;i<num_dsegment;i++) // the flow is simplified as one-dimension flow
		{
			SgtPropUpdate(segpoint, Dehumidifier);
			double fa_rho_value=segpoint.fa_rho;
			double fa_cp_value=segpoint.fa_cp;	
			double fa_htc_value=segpoint.fa_htc;
			double fa_mtc_value=segpoint.fa_mtc;
			double fa_tin_value=segpoint.fa_temp;
			double fa_flin_value = segpoint.fa_flow;
			double fa_vapin_value=segpoint.fa_vap;
			double pf_rho_value=segpoint.pf_rho;
			double pf_cp_value=segpoint.pf_cp;	
			double pf_htc_value=segpoint.pf_htc;
			double pf_mtc_value=segpoint.pf_mtc;
			double pf_tin_value=segpoint.pf_temp;
			double pf_flin_value=segpoint.pf_flow;
			double pf_vapin_value=segpoint.pf_vap;
			double kt_mem_value=segpoint.m_kt;

			double pvap_fa = vapPstaturate(temp_fa) * airrelhum(temp_fa, humid_fa);
			double pvap_pf=0.0;
			double mfl_mem_value = 0.0;
			if(proppf_airvap == 0) {pvap_pf = press_pf;}
			if(proppf_airvap == 1) {pvap_pf = vapPstaturate(temp_pf) * airrelhum(temp_pf, humid_pf);}
			if(prop_memmodel == 0)
			{
				double pp_mvap = (pvap_fa + pvap_pf) / 2.0;
				double pp_mair = segpoint.fa_p - pvap_fa;
				mfl_mem_value = watpermeatd(segpoint, Dehumidifier, pvap_fa, pvap_pf, temp_fa, pp_mvap, pp_mair);
			}
			if (prop_memmodel == 1)
			{
				double rh_fa=airrelhum(temp_fa, humid_fa);
				double t_fa=temp_fa;
				double fl_fa=fa_flin_value;
				double d_fa = Dehumidifier->get_lcseg_fa()/(Dehumidifier->aseg_crossfa());
				mfl_mem_value = watpermeatdais(pvap_fa, pvap_pf, rh_fa, fl_fa, t_fa, d_fa);
			}
			if (prop_memmodel == 2)
			{mfl_mem_value = watpermeatc(pvap_fa, pvap_pf);}
			
			SegmentSolve sgt(fa_flin_value, fa_rho_value, fa_cp_value, lseg_fa_value, dseg_fa_value, fa_htc_value, fa_mtc_value, fa_tin_value, fa_vapin_value,
							 pf_flin_value, pf_rho_value, pf_cp_value, lseg_pf_value, dseg_pf_value, pf_htc_value, pf_mtc_value, pf_tin_value, pf_vapin_value,
							 delta_mem_value, kt_mem_value, mfl_mem_value, ftype, laynvalue, enable_mtimpact,memsurf_correct);

			temp_fa=sgt.get_fa_temp();
			flow_fa=flow_fa-mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_fa_coeff()*memsurf_correct;
			vap_fa=sgt.get_fa_vap();
			humid_fa=vap_fa/(mflwseg_fa_value-mflwvapseg_fa_value);
			press_fa=segpoint.fa_p; 
			temp_pf=sgt.get_pf_temp();
			flow_pf=flow_pf+mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_pf_coeff()*memsurf_correct;
			vap_pf=sgt.get_pf_vap();
			if(proppf_airvap == 0) {humid_pf=1.0;}
			if (proppf_airvap == 0) { temp_pf = sgt.get_pfm_temp(); }
			if(proppf_airvap == 1) {humid_pf=vap_pf/(mflwseg_pf_value-mflwvapseg_pf_value);}
			press_pf=segpoint.pf_p;

			map_temp_fa[0][i]=temp_fa;
			map_temp_fam[0][i]=sgt.get_fam_temp();
			map_vap_fa[0][i]=vap_fa;
			map_humid_fa[0][i]= humid_fa;
			map_pvap_fa[0][i]=pvap_fa;
			map_flow_fa[0][i]=flow_fa;
			map_temp_pf[0][i]=temp_pf;
			map_temp_pfm[0][i]=sgt.get_pfm_temp();
			map_vap_pf[0][i]=vap_pf;
			map_humid_pf[0][i]= humid_pf;
			map_pvap_pf[0][i]=pvap_pf;
			map_flow_pf[0][i]=flow_pf;
			map_permeat_vap[0][i]= mfl_mem_value*memsurf_correct;
		}
		tfa_out= map_temp_fa[0][num_dsegment-1];
		vapfa_out= map_vap_fa[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_fa));
		flfa_out= map_flow_fa[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_fa));
		tpf_out= map_temp_pf[0][num_dsegment - 1];
		vappf_out= map_vap_pf[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_pf));
		flpf_out= map_flow_pf[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_pf));
	}
	else if (flowtype == 2)
	{
		//             fa_in^pf_out      : parallel counterflow
		//			 ______|____|____
		//		Hpf	|      |    |   |\ 
		//		Hm	|======v========| \
		//		Hfa	|_______________|\ \
		//			\      L        \D\|
		//			 \           ^   \ |
		//			  \__________|____\|
		//				   |     | pf_in
		//				   |
		//				   V fa_out
		//
		//cout<<"Notice: it is a parallel counterflow along D dirction"<<endl;
		double pf_coeff=0.0;
		if(num_layer==1){pf_coeff=1.0;}
		else if(num_layer==2){pf_coeff=2.0;}
		else if(num_layer>2){pf_coeff=2.0;}

		double terror=0.0;
		double verror=0.0;
		double tmin=min(tfa_in,tpf_in)-0.01;
		double tmax=max(tfa_in,tpf_in)+0.01;
		double vmin=min(mflwvapseg_fa_value,mflwvapseg_pf_value)*pf_coeff*0.0;
		double vmax=max(mflwvapseg_fa_value,mflwvapseg_pf_value)*pf_coeff*1.0;
		double tcmd=0.0;
		double vcmd=0.0;
		double cumvfl=0.0;
		double omvfl=0.0;
		int ido=0;
		do{
			ido++;
			tcmd=(tmin+tmax)/2.0;
			vcmd=(vmin+vmax)/2.0;

			temp_fa = tfa_in;
			flow_fa = mflwseg_fa_value;
			vap_fa= mflwvapseg_fa_value;
			humid_fa=mflwvapseg_fa_value/(mflwseg_fa_value-mflwvapseg_fa_value);
			press_fa= pfa_in; 
			temp_pf = tcmd;
			flow_pf = mflwseg_pf_value+vcmd;
			vap_pf = mflwvapseg_pf_value+vcmd;
			if(proppf_airvap == 0) {humid_pf=1.0;}
			if(proppf_airvap == 1) {humid_pf=vap_pf/(mflwseg_pf_value-mflwvapseg_pf_value);}
			press_pf= ppf_in;
			cumvfl	= 0.0;

			for(int i=0;i<num_dsegment;i++) // the flow is simplified as one-dimension flow
			{
				SgtPropUpdate(segpoint, Dehumidifier);
				double fa_rho_value=segpoint.fa_rho;
				double fa_cp_value=segpoint.fa_cp;
				double fa_htc_value=segpoint.fa_htc;
				double fa_mtc_value=segpoint.fa_mtc;
				double fa_tin_value=segpoint.fa_temp;
				double fa_flin_value=segpoint.fa_flow;
				double fa_vapin_value=segpoint.fa_vap;
				double pf_rho_value=segpoint.pf_rho;
				double pf_cp_value=segpoint.pf_cp;
				double pf_htc_value=segpoint.pf_htc;
				double pf_mtc_value=segpoint.pf_mtc;
				double pf_tin_value=segpoint.pf_temp;
				double pf_flin_value=segpoint.pf_flow;
				double pf_vapin_value=segpoint.pf_vap;
				double kt_mem_value=segpoint.m_kt;
				double pvap_fa = vapPstaturate(temp_fa) * airrelhum(temp_fa, humid_fa);
				double pvap_pf =0;
				double mfl_mem_value =0;
				if(proppf_airvap == 0) {pvap_pf = press_pf;}
				if(proppf_airvap == 1) {pvap_pf = vapPstaturate(temp_pf) * airrelhum(temp_pf, humid_pf);}
				if (prop_memmodel == 0)
				{
					double pp_mvap = (pvap_fa + pvap_pf) / 2.0;
					double pp_mair = segpoint.fa_p - pvap_fa;
					mfl_mem_value = watpermeatd(segpoint, Dehumidifier, pvap_fa, pvap_pf, temp_fa, pp_mvap, pp_mair);
				}
				if (prop_memmodel == 1)
				{
					double rh_fa = airrelhum(temp_fa, humid_fa);
					double t_fa = temp_fa;
					double fl_fa = fa_flin_value;
					double d_fa = Dehumidifier->get_lcseg_fa() / (Dehumidifier->aseg_crossfa());
					mfl_mem_value = watpermeatdais(pvap_fa, pvap_pf, rh_fa, fl_fa, t_fa, d_fa);
				}
				if (prop_memmodel == 2)
				{mfl_mem_value = watpermeatc(pvap_fa, pvap_pf);}

				SegmentSolve sgt(fa_flin_value, fa_rho_value, fa_cp_value, lseg_fa_value, dseg_fa_value, fa_htc_value, fa_mtc_value, fa_tin_value, fa_vapin_value,
							     pf_flin_value, pf_rho_value, pf_cp_value, lseg_pf_value, dseg_pf_value, pf_htc_value, pf_mtc_value, pf_tin_value, pf_vapin_value,
								 delta_mem_value, kt_mem_value, mfl_mem_value, ftype, laynvalue, enable_mtimpact,memsurf_correct);

				cumvfl=cumvfl+mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_pf_coeff()*memsurf_correct;
				temp_fa=sgt.get_fa_temp();
				flow_fa=flow_fa-mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_fa_coeff()*memsurf_correct;
				vap_fa=sgt.get_fa_vap();
				humid_fa=vap_fa/(mflwseg_fa_value-mflwvapseg_fa_value);
				press_fa=segpoint.fa_p;
				temp_pf=sgt.get_pf_temp();
				flow_pf=flow_pf-mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_pf_coeff()*memsurf_correct;
				if (i == 0) {omvfl=mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_pf_coeff()*memsurf_correct;}
				vap_pf=sgt.get_pf_vap();
				if(proppf_airvap == 0) {humid_pf=1.0;temp_pf=sgt.get_pfm_temp();}
				if(proppf_airvap == 1) {humid_pf=vap_pf/(mflwseg_pf_value-mflwvapseg_pf_value);}
				press_pf=segpoint.pf_p;

				map_temp_fa[0][i]=temp_fa;
				map_temp_fam[0][i]=sgt.get_fam_temp();
				map_vap_fa[0][i]=vap_fa;
				map_humid_fa[0][i]=humid_fa;
				map_pvap_fa[0][i]=pvap_fa;
				map_flow_fa[0][i]=flow_fa;
				map_temp_pf[0][i]=temp_pf;
				map_temp_pfm[0][i]= sgt.get_pfm_temp();
				map_vap_pf[0][i]=vap_pf;
				map_humid_pf[0][i]=humid_pf;
				map_pvap_pf[0][i]=pvap_pf;
				map_flow_pf[0][i]=flow_pf;
				map_permeat_vap[0][i] = mfl_mem_value*memsurf_correct;
			}
			terror=(tpf_in-map_temp_pf[0][num_dsegment-1]);
			if (abs(terror)>0.01)
			{
				if (terror>0.01)
				{tmin = tcmd;}
				else if (terror<-0.01)
				{tmax = tcmd;}
			}

			verror=(1-cumvfl/vcmd);
			if (abs(verror)>0.01)
			{
				if(verror>=0.01)
				{vmax=vcmd;}
				else if(verror<=-0.01)
				{vmin=vcmd;}
			}
			//cout << "Loop running no. " << ido << "; Terror:" << terror << "; Verror:" << verror <<endl;
			if (ido > 50 && abs(terror)>0.01) {tmin=tmin+terror/10.0; tmax=tmax+terror/10.0;}
			if (ido > 50 && abs(verror)>0.01) {vmin=vmin+verror/10.0; vmax=vmax+verror/10.0;}
			if (ido > 500) break;
		}while(abs(verror)>0.01 || abs(terror)>0.01);
		tfa_out = map_temp_fa[0][num_dsegment - 1];
		vapfa_out = map_vap_fa[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_fa));
		flfa_out = map_flow_fa[0][num_dsegment - 1]*(static_cast<double>(num_lsegment*numchannel_fa));
		tpf_out = map_temp_pf[0][0];
		vappf_out = (map_vap_pf[0][0]+ omvfl)*(static_cast<double>(num_lsegment*numchannel_pf));
		flpf_out = (map_flow_pf[0][0]+ omvfl)*(static_cast<double>(num_lsegment*numchannel_pf));
	}
	else
	{
		//             fa_int			: crossflow
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
		//	
		//cout<<"Notice: it is a crossflow"<<endl; 
		temp_fa=tfa_in;
		flow_fa= mflwseg_fa_value;
		vap_fa= mflwvapseg_fa_value;
		humid_fa=mflwvapseg_fa_value/(mflwseg_fa_value-mflwvapseg_fa_value);
		press_fa=pfa_in; 
		temp_pf=tpf_in;
		flow_pf = mflwseg_pf_value;
		vap_pf=mflwvapseg_pf_value;
		if(proppf_airvap == 0) {humid_pf=1.0;}
		if(proppf_airvap == 1) {humid_pf=mflwvapseg_pf_value/(mflwseg_pf_value-mflwvapseg_pf_value);}
		press_pf=ppf_in;

		for(int i=0;i<num_lsegment;i++)
		{
				for(int j=0;j<num_dsegment;j++)
				{
					if(i<1)
					{
						temp_pf=tpf_in;
						flow_pf= mflwseg_pf_value;
						vap_pf=mflwvapseg_pf_value;
						if(proppf_airvap == 0) {humid_pf=1.0;}
						if(proppf_airvap == 1) {humid_pf=mflwvapseg_pf_value/(mflwseg_pf_value-mflwvapseg_pf_value);}
						press_pf=ppf_in;
					}
					else
					{
						temp_pf=map_temp_pf[i-1][j];
						flow_pf=map_flow_pf[i-1][j];
						vap_pf=map_vap_pf[i-1][j];
						if(proppf_airvap == 0) {humid_pf=1.0;}
						if(proppf_airvap == 1) {humid_pf=vap_pf/(mflwseg_pf_value-mflwvapseg_pf_value);}
						press_pf=ppf_in;

					}
					if(j<1)
					{
						temp_fa=tfa_in;
						flow_fa= mflwseg_fa_value;
						vap_fa=mflwvapseg_fa_value;
						humid_fa=mflwvapseg_fa_value/(mflwseg_fa_value-mflwvapseg_fa_value);
						press_fa=pfa_in;
					}
					SgtPropUpdate(segpoint, Dehumidifier);
					double fa_rho_value=segpoint.fa_rho;
					double fa_cp_value=segpoint.fa_cp;
					double fa_htc_value=segpoint.fa_htc;
					double fa_mtc_value=segpoint.fa_mtc;
					double fa_tin_value=segpoint.fa_temp;
					double fa_flin_value=segpoint.fa_flow;
					double fa_vapin_value=segpoint.fa_vap;
					double pf_rho_value=segpoint.pf_rho;
					double pf_cp_value=segpoint.pf_cp;
					double pf_htc_value=segpoint.pf_htc;
					double pf_mtc_value=segpoint.pf_mtc;
					double pf_tin_value=segpoint.pf_temp;
					double pf_flin_value=segpoint.pf_flow;
					double pf_vapin_value=segpoint.pf_vap;
					double kt_mem_value=segpoint.m_kt;

					double pvap_fa =vapPstaturate(temp_fa)*airrelhum(temp_fa, humid_fa);
					double pvap_pf =0.0;
					double mfl_mem_value = 0.0;
					if(proppf_airvap == 0) {pvap_pf = press_pf;}
					if(proppf_airvap == 1) {pvap_pf = vapPstaturate(temp_pf) * airrelhum(temp_pf, humid_pf);}
					if (prop_memmodel == 0)
					{
						double pp_mvap = (pvap_fa + pvap_pf) / 2.0;
						double pp_mair = segpoint.fa_p - pvap_fa;
						mfl_mem_value = watpermeatd(segpoint, Dehumidifier, pvap_fa, pvap_pf, temp_fa, pp_mvap, pp_mair);
					}
					if (prop_memmodel == 1)
					{
						double rh_fa = airrelhum(temp_fa, humid_fa);
						double t_fa = temp_fa;
						double fl_fa = fa_flin_value;
						double d_fa = Dehumidifier->get_lcseg_fa() / (Dehumidifier->aseg_crossfa());
						mfl_mem_value = watpermeatdais(pvap_fa, pvap_pf, rh_fa, fl_fa, t_fa, d_fa);
					}
					if (prop_memmodel == 2)
					{mfl_mem_value = watpermeatc(pvap_fa, pvap_pf);}

					SegmentSolve sgt(fa_flin_value, fa_rho_value, fa_cp_value, lseg_fa_value, dseg_fa_value, fa_htc_value, fa_mtc_value, fa_tin_value, fa_vapin_value,
									 pf_flin_value, pf_rho_value, pf_cp_value, lseg_pf_value, dseg_pf_value, pf_htc_value, pf_mtc_value, pf_tin_value, pf_vapin_value,
									 delta_mem_value, kt_mem_value, mfl_mem_value, ftype, laynvalue, enable_mtimpact,memsurf_correct);

					temp_fa=sgt.get_fa_temp();
					flow_fa=flow_fa-mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_fa_coeff()*memsurf_correct;
					vap_fa=sgt.get_fa_vap();
					humid_fa=vap_fa/(mflwseg_fa_value-mflwvapseg_fa_value);
					press_fa=segpoint.fa_p;
					temp_pf=sgt.get_pf_temp();
					flow_pf=flow_pf+mfl_mem_value*dseg_fa_value*lseg_fa_value*sgt.ask_pf_coeff()*memsurf_correct;
					vap_pf=sgt.get_pf_vap();
					if(proppf_airvap == 0) {humid_pf=1.0;}
					if(proppf_airvap == 1) {humid_pf=vap_pf/(mflwseg_pf_value-mflwvapseg_pf_value);}
					press_pf=segpoint.pf_p;

					map_temp_fa[i][j]=temp_fa;
					map_temp_fam[i][j]=sgt.get_fam_temp();
					map_vap_fa[i][j]=vap_fa;
					map_humid_fa[i][j]=humid_fa;
					map_pvap_fa[i][j]= pvap_fa;
					map_flow_fa[i][j]=flow_fa;
					map_temp_pf[i][j]=temp_pf;
					map_temp_pfm[i][j]=sgt.get_pfm_temp();
					map_vap_pf[i][j]=vap_pf;
					map_humid_pf[i][j]=humid_pf;
					map_pvap_pf[i][j]= pvap_pf;
					map_flow_pf[i][j]=flow_pf;
					map_permeat_vap[i][j] = mfl_mem_value*memsurf_correct;
				}
		}
		tfa_out = 0;
		vapfa_out = 0;
		flfa_out = 0;
		tpf_out = 0;
		vappf_out = 0;
		flpf_out = 0;
		for (int k = 0; k < num_lsegment; k++)
		{
			tfa_out = tfa_out+map_temp_fa[k][num_dsegment-1]/(static_cast<double>(num_lsegment));
			vapfa_out = vapfa_out+map_vap_fa[k][num_dsegment - 1]*(static_cast<double>(numchannel_fa));
			flfa_out = flfa_out+map_flow_fa[k][num_dsegment - 1]*(static_cast<double>(numchannel_fa));
		}
		for (int k = 0; k < num_dsegment; k++)
		{
			tpf_out = tpf_out+map_temp_pf[num_lsegment-1][k]/(static_cast<double>(num_dsegment));
			vappf_out = vappf_out+map_vap_pf[num_lsegment-1][k]*(static_cast<double>(numchannel_pf));
			flpf_out = flpf_out+map_flow_pf[num_lsegment-1][k]*(static_cast<double>(numchannel_pf));
		}
	}
}//end of energy/mass analysis for all segments of a component device

void Component::SgtPropUpdate(SegmentProperty &segpoint, HXM *hExchanger)
{
		segpoint.m_rp=rp_membrane;
		segpoint.m_ps=ps_membrane;
		segpoint.m_ts=ts_membrane;
		segpoint.m_kt=kt_membrane;
		segpoint.fa_temp=temp_fa;
		segpoint.fa_flow=flow_fa;
		segpoint.fa_vap=vap_fa;
		segpoint.fa_hr=humid_fa;
		segpoint.fa_p=press_fa;
		segpoint.pf_temp=temp_pf;
		segpoint.pf_flow=flow_pf;
		segpoint.pf_vap=vap_pf;
		segpoint.pf_hr=humid_pf;
		segpoint.pf_p=press_pf;
		segpoint.prop_airvap = proppf_airvap;
		segpoint.htc_correct = htc_correction;
		segpoint.ftc_correct = memsurf_correction;
		segpoint.memsurf_correct = ftc_correction;
		allprop_flatplate(segpoint, hExchanger, 1);
		allprop_flatplate(segpoint, hExchanger, 2);
}//end of udpating segmentsProperty
