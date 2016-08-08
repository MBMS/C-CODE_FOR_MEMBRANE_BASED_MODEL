//SegmentSolve.cpp written by Zhiming Gao (gaoz@ornl.gov)


#include <iostream>
#include <cstdlib>
#include "SegmentSolve.h"

using namespace std;
using namespace MoistAirVapProperty;


//This is a function that create the matrix A and call jfc for solving HT in the segment written by Zhiming Gao
SegmentSolve::SegmentSolve(double mflwseg_fa_value, double fa_rho_value, double fa_cp_value, double lseg_fa_value, double dseg_fa_value, 
						   double fa_htc_value, double fa_mtc_value, double fa_tin_value, double fa_vapin_value,
						   double mflwseg_pf_value, double pf_rho_value, double pf_cp_value, double lseg_pf_value, double dseg_pf_value, 
						   double pf_htc_value, double pf_mtc_value, double pf_tin_value, double pf_vapin_value,
						   double delta_mem_value, double kt_mem_value, double mfl_mem_value, int ftype, int mem_laynvalue, 
						   double enable_mtimpact, double memsurf_correct)
{
	mflwseg_fa=mflwseg_fa_value;	
	fa_rho=fa_rho_value;	
	fa_cp=fa_cp_value;	
	lseg_fa=lseg_fa_value;
	dseg_fa=dseg_fa_value;
	fa_htc=fa_htc_value;	
	fa_mtc=fa_mtc_value;
	fa_tin=fa_tin_value;	
	fa_vapin=fa_vapin_value;	
	mflwseg_pf=mflwseg_pf_value;
	pf_rho=pf_rho_value;	
	pf_cp=pf_cp_value;	
	lseg_pf=lseg_pf_value;
	dseg_pf=dseg_pf_value;
	pf_htc=pf_htc_value;	
	pf_mtc=pf_mtc_value;
	pf_tin=pf_tin_value;	
	pf_vapin=pf_vapin_value;		
	delta_mem=delta_mem_value;	
	kt_mem=kt_mem_value;		
	mfl_mem=mfl_mem_value;	
	flowtype=ftype;
	enable_mtm = enable_mtimpact;
	enable_msf= memsurf_correct;
	num_membrane=mem_laynvalue;

	get_fa_coeff();
	get_pf_coeff();

	SegmentSolht();
	SegmentSolmt();
}


void SegmentSolve::SegmentSolht()
{
	double eps=1.0e-10;
	double a[4][4]={0.0};
	double b[4] = {0.0};
	int varsize=sizeof(b)/sizeof(b[0]);

	double dx = dseg_fa;
	double dy = lseg_fa;
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
		double mcpdx_fa=mflwseg_fa*fa_cp/dx;
		double htcdy2_fa=fa_coeff*fa_htc*dy*enable_msf;

		double mcpdx_pf=mflwseg_pf*pf_cp/dx;
		double htcdy2_pf=pf_coeff*pf_htc*dy*enable_msf;

		double htc_fa=fa_htc;
		double htc_pf=pf_htc;
		double ktdelta_mem=kt_mem/delta_mem;

		double fwcp_fa = fa_coeff*mfl_mem*dy*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;
		double fwcp_fam= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pfm= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pf = pf_coeff*mfl_mem*dy*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;

		a[0][0]=mcpdx_fa+htcdy2_fa/2.0+fwcp_fa/2.0; a[0][1]=(-1.0)*htcdy2_fa-fwcp_fa; b[0]=(mcpdx_fa-htcdy2_fa/2.0-fwcp_fa/2.0)*fa_tin;
		a[1][0]=htc_fa/2.0+fwcp_fam/2.0; a[1][1]=(-1.0)*(ktdelta_mem+htc_fa); a[1][2]=ktdelta_mem-fwcp_fam; b[1]=(-1.0)*(htc_fa/2.0+fwcp_fam/2.0)*fa_tin;
		a[2][1]=(-1.0)*(ktdelta_mem+fwcp_pfm); a[2][2]=(ktdelta_mem+htc_pf); a[2][3]=(-1.0)*htc_pf/2.0+fwcp_pfm/2.0; b[2]=(htc_pf/2.0-fwcp_pfm/2.0)*pf_tin;
		a[3][2]=(-1.0)*(htcdy2_pf)-fwcp_pf; a[3][3]=mcpdx_pf+htcdy2_pf/2.0+fwcp_pf/2.0; b[3]=(mcpdx_pf-htcdy2_pf/2.0-fwcp_pf/2.0)*pf_tin;

	}
	else if (flowtype == 2)
	{
		//            fa_in^pf_out      : parallel counterflow
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
		double mcpdx_fa=mflwseg_fa*fa_cp/dx;
		double htcdy2_fa=fa_coeff*fa_htc*dy*enable_msf;

		double mcpdx_pf=mflwseg_pf*pf_cp/dx;
		double htcdy2_pf=pf_coeff*pf_htc*dy*enable_msf;

		double htc_fa=fa_htc;
		double htc_pf=pf_htc;
		double ktdelta_mem=kt_mem/delta_mem;

		double fwcp_fa = fa_coeff*mfl_mem*dy*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;
		double fwcp_fam= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pfm= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pf = pf_coeff*mfl_mem*dy*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;

		a[0][0]=mcpdx_fa+htcdy2_fa/2.0+fwcp_fa/2.0; a[0][1]=(-1.0)*htcdy2_fa-fwcp_fa; b[0]=(mcpdx_fa-htcdy2_fa/2.0-fwcp_fa/2.0)*fa_tin;
		a[1][0]=htc_fa/2.0+fwcp_fam/2.0; a[1][1]=(-1.0)*(ktdelta_mem+htc_fa); a[1][2]=ktdelta_mem-fwcp_fam; b[1]=(-1.0)*(htc_fa/2.0+fwcp_fam/2.0)*fa_tin;
		a[2][1]=(-1.0)*(ktdelta_mem+fwcp_pfm); a[2][2]=(ktdelta_mem+htc_pf); a[2][3]=(-1.0)*htc_pf/2.0+fwcp_pfm/2.0; b[2]=(htc_pf/2.0-fwcp_pfm/2.0)*pf_tin;
		a[3][2]=htcdy2_pf+fwcp_pf; a[3][3]=mcpdx_pf-htcdy2_pf/2.0-fwcp_pf/2.0; b[3]=(mcpdx_pf+htcdy2_pf/2.0+fwcp_pf/2.0)*pf_tin;

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
		double mcpdx_fa=mflwseg_fa*fa_cp/dx;
		double htcdy2_fa=fa_coeff*fa_htc*dy*enable_msf;

		double mcpdy_pf=mflwseg_pf*pf_cp/dy;
		double htcdx2_pf=pf_coeff*pf_htc*dx*enable_msf;

		double htc_fa=fa_htc;
		double htc_pf=pf_htc;
		double ktdelta_mem=kt_mem/delta_mem;

		double fwcp_fa = fa_coeff*mfl_mem*dy*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;
		double fwcp_fam= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pfm= mfl_mem*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm;
		double fwcp_pf = pf_coeff*mfl_mem*dx*vapcp((fa_tin+pf_tin)/2.0)*enable_mtm*enable_msf;

		a[0][0]=mcpdx_fa+htcdy2_fa/2.0+fwcp_fa/2.0; a[0][1]=(-1.0)*htcdy2_fa-fwcp_fa; b[0]=(mcpdx_fa-htcdy2_fa/2.0-fwcp_fa/2.0)*fa_tin;
		a[1][0]=htc_fa/2.0+fwcp_fam/2.0; a[1][1]=(-1.0)*(ktdelta_mem+htc_fa); a[1][2]=ktdelta_mem-fwcp_fam; b[1]=(-1.0)*(htc_fa/2.0+fwcp_fam/2.0)*fa_tin;
		a[2][1]=(-1.0)*(ktdelta_mem+fwcp_pfm); a[2][2]=(ktdelta_mem+htc_pf); a[2][3]=(-1.0)*htc_pf/2.0+fwcp_pfm/2.0; b[2]=(htc_pf/2.0-fwcp_pfm/2.0)*pf_tin;
		a[3][2]=(-1.0)*htcdx2_pf-fwcp_pf; a[3][3]=mcpdy_pf+htcdx2_pf/2.0+fwcp_pf/2.0; b[3]=(mcpdy_pf-htcdx2_pf/2.0-fwcp_pf/2.0)*pf_tin;

	}

	jfc(a, b, varsize, eps);
	fa_temp_value=b[0];
	fam_temp_value=b[1];
	pfm_temp_value=b[2];
	pf_temp_value=b[3]; //for parallel :pf_temp_value is segment-out; for parallel counterflow: pf_temp_value is segment-in; for crossflow: segment-out;
} //end of SegmentSolve for energy equations


void SegmentSolve::SegmentSolmt()
{
	double eps=1.0e-10;
	double fa_vapout=0.0;
	double fam_vap=0.0;
	double pf_vapout=0.0;
	double pfm_vap=0.0;

	double dx = dseg_fa;
	double dy = lseg_fa;
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
		double mtdydx_fa=fa_coeff*mfl_mem*dy*dx*enable_msf;
		fa_vapout=fa_vapin-mtdydx_fa;
		fam_vap=mtdydx_fa;

		double mtdydx_pf=pf_coeff*mfl_mem*dy*dx*enable_msf;
		pf_vapout=pf_vapin+mtdydx_pf; 
		pfm_vap=mtdydx_fa;
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
		double mtdydx_fa=fa_coeff*mfl_mem*dy*dx*enable_msf;
		fa_vapout=fa_vapin-mtdydx_fa;
		fam_vap=mtdydx_fa;

		double mtdydx_pf=pf_coeff*mfl_mem*dy*dx*enable_msf;
		pf_vapout=pf_vapin-mtdydx_pf; 
		pfm_vap=mtdydx_fa;
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
		double mtdydx_fa=fa_coeff*mfl_mem*dy*dx*enable_msf;
		fa_vapout=fa_vapin-mtdydx_fa;
		fam_vap=mtdydx_fa;

		double mtdydx_pf=pf_coeff*mfl_mem*dy*dx*enable_msf;
		pf_vapout=pf_vapin+mtdydx_pf; 
		pfm_vap=mtdydx_fa;
	}
	fa_vap_value=fa_vapout;
	fam_vap_value=fam_vap;
	pfm_vap_value=pfm_vap;
	pf_vap_value=pf_vapout; //for parallel :pf_w_value is segment-out; for parallel counterflow: is segment-in; for crossflow: segment-out;
	if(fa_vap_value<0) fa_vap_value=0;
	if(fam_vap_value<0) fam_vap_value=0;
	if(pf_vap_value<0) pf_vap_value=0;
	if(pfm_vap_value<0) pfm_vap_value=0;
} //end of SegmentSolms for mass equations

//This is a function that calculate the matrix A*X=B written by Zhiming Gao
void SegmentSolve::jfc(double a[4][4], double *b, int m, double eps)
{	
	double c=0.0;
	double cmax=0.0;
	double cmax_org=0.0; 
	double cmd =0.0;

	for(int i=0;i<m-1;i++)
	{
		cmax=abs(a[i][i]);
		cmax_org=a[i][i];
		for(int j=i+1;j<m;j++)
		{
			if(abs(a[j][i])>cmax)
			{
				cmax= abs(a[j][i]);
				cmax_org=a[j][i];
				for(int k=i;k<m;k++)
				{
				cmd=a[i][k];
				a[i][k]=a[j][k];
				a[j][k]=cmd;
				}
				cmd=b[i];
				b[i]=b[j];
				b[j]=cmd;
			}
		}

		if(cmax<eps){cout<<"warming: max<eps"<<endl;}
		for(int j1=i;j1<m;j1++)
		{a[i][j1]=a[i][j1]/cmax_org;}
		b[i]=b[i]/cmax_org;

		for(int k1=i+1;k1<m;k1++)
		{
			c=a[k1][i];
			for(int n1=i;n1<m;n1++)
			{a[k1][n1]=a[k1][n1]-a[i][n1]*c;}
			b[k1]=b[k1]-b[i]*c;
		}
	}
	b[m-1]=b[m-1]/a[m-1][m-1];

	for(int j2=m-1;j2>-1;j2--)
	{
		for(int k2=j2-1;k2>-1;k2--)
		{b[k2]=b[k2]-b[j2]*a[k2][j2];}
	}
} //end of Gaussian Elimination


double SegmentSolve::get_fa_temp()	//get fa_temp
{return fa_temp_value;}

double SegmentSolve::get_fam_temp()	//get fam_temp
{return fam_temp_value;}

double SegmentSolve::get_pfm_temp()	//get pfm_temp
{return pfm_temp_value;}

double SegmentSolve::get_pf_temp()	//get pf_temp
{return pf_temp_value;}

double SegmentSolve::get_fa_vap()	//get fa_w
{return fa_vap_value;}

double SegmentSolve::get_fam_vap()	//get fam_w
{return fam_vap_value;}

double SegmentSolve::get_pfm_vap()	//get pfm_w
{return pfm_vap_value;}

double SegmentSolve::get_pf_vap()	//get pf_w
{return pf_vap_value;}

void SegmentSolve::get_fa_coeff()	//get fa_coeff
{
	if(num_membrane ==1){fa_coeff=1.0;}
	else if(num_membrane ==2){fa_coeff=1.0;}
	else if(num_membrane >2){fa_coeff=2.0;}
	else {cout<<"The value for membrane layer is wrong, please check!"<<endl;system("pause");exit(1);}
	//return fa_coeff;
}

void SegmentSolve::get_pf_coeff()	//get pf_coeff
{
	if(num_membrane ==1){pf_coeff=1.0;}
	else if(num_membrane ==2){pf_coeff=2.0;}
	else if(num_membrane >2){pf_coeff=2.0;}
	else {cout<<"The value for membrane layer is wrong, please check!"<<endl;system("pause");exit(1);}
	//return pf_coeff;
}

double SegmentSolve::ask_fa_coeff()	//get fa_coeff
{return fa_coeff;}

double SegmentSolve::ask_pf_coeff()	//get pf_coeff
{return pf_coeff;}