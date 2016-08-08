//HXM.cpp written by Zhiming Gao (gaoz@ornl.gov)


#include <iostream>
#include <cstdlib>
#include "HXM.h"

using namespace std;


HXM::HXM(double lvalue, double dvalue, double tfavalue, double tpfvalue, double tmvalue, double mffa, double mfvapfa, double mfpf, double mfvappf, double memsurfr, int lnvalue, int dnvalue, int laynvalue, int ftype)
{
		length=lvalue;
		depth =dvalue;
		thick_fa=tfavalue; 
		thick_pf=tpfvalue;
		thick_membrane=tmvalue;
		mf_fa=mffa;
		mfvap_fa=mfvapfa;
		mf_pf=mfpf;
		mfvap_pf=mfvappf;
		memsur_corr=memsurfr;
		num_lsegment=lnvalue;
		num_dsegment=dnvalue;
		num_layer=laynvalue;
		flowtype=ftype;			// 1: parallel uncounterflow along D dirction; 2:parallel counterflow along D dirction; 0 or other: crossflow

		numchannel_fa=num_layer/2+1;
		numchannel_pf=(num_layer+1)/2;
		if (num_layer>2) {numchannel_fa=(num_layer+1)/2; numchannel_pf=(num_layer+1)/2;}
		if (flowtype == 1)
		{
			/*             fa_in|pf_in      : parallel uncounterflow
					 ______|____|____
				Hpf	|      |    v   |\ 
				Hm	|======v========| \
				Hfa	|_______________|\ \
					\      L        \D\|
					 \           |   \ |
					  \__________|____\|
						   |     V pf_out
						   |
						   V fa_out
			*/ 
			cout<<"Notice: it is a parallel uncounterflow along D dirction"<<endl; 
			lseg_fa=length/(static_cast<double>(num_lsegment));	
			dseg_fa=depth/(static_cast<double>(num_dsegment));	
			mflwseg_fa=mf_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			mflwvapseg_fa=mfvap_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			lseg_pf=lseg_fa; 
			dseg_pf=dseg_fa;	
			mflwseg_pf=mf_pf/(static_cast<double>(num_lsegment*numchannel_pf)); 
			mflwvapseg_pf=mfvap_pf/(static_cast<double>(num_lsegment*numchannel_pf)); 
		}
		else if (flowtype == 2)
		{
			/*             fa_in^pf_out      : parallel counterflow
					 ______|____|____
				Hpf	|      |    |   |\ 
				Hm	|======v========| \
				Hfa	|_______________|\ \
					\      L        \D\|
					 \           ^   \ |
					  \__________|____\|
						   |     | pf_in
						   |
						   V fa_out
			*/
			cout<<"Notice: it is a parallel counterflow along D dirction"<<endl; 
			lseg_fa=length/(static_cast<double>(num_lsegment));	
			dseg_fa=depth/(static_cast<double>(num_dsegment));	
			mflwseg_fa=mf_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			mflwvapseg_fa=mfvap_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			lseg_pf=lseg_fa; 
			dseg_pf=dseg_fa;	
			mflwseg_pf=mf_pf/(static_cast<double>(num_lsegment*numchannel_pf)); 
			mflwvapseg_pf=mfvap_pf/(static_cast<double>(num_lsegment*numchannel_pf)); 
		}
		else
		{
			/*             fa_int			: crossflow
					 ______|________
			-->	Hpf	|      |        |\ 
				Hm	|======v========| \
				Hfa	|_______________|\-\-> pf_out
					\      L        \D\|
					 \               \ |
					  \_______________\|
							   |
							   |
							   V fa_out
			*/
			cout<<"Notice: it is a crossflow"<<endl; 
			lseg_fa=length/(static_cast<double>(num_lsegment));	
			dseg_fa=depth/(static_cast<double>(num_dsegment));	
			mflwseg_fa=mf_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			mflwvapseg_fa=mfvap_fa/(static_cast<double>(num_lsegment*numchannel_fa));	
			lseg_pf=dseg_fa; 
			dseg_pf=lseg_fa;	
			mflwseg_pf=mf_pf/(static_cast<double>(num_dsegment*numchannel_pf)); 
			mflwvapseg_pf=mfvap_pf/(static_cast<double>(num_dsegment*numchannel_pf)); 
		}
}

double HXM::aseg_ht()
{return lseg_fa*dseg_fa;}

double HXM::aseg_crossfa()
{return lseg_fa*thick_fa;}

double HXM::aseg_crosspf()
{return lseg_pf*thick_pf;}

double HXM::get_lseg_fa()
{return lseg_fa;}

double HXM::get_dseg_fa()
{return dseg_fa;}

double HXM::get_mflwseg_fa()
{return mflwseg_fa;}

double HXM::get_mflwvapseg_fa()
{return mflwvapseg_fa;}

double HXM::get_lseg_pf()
{return lseg_pf;}

double HXM::get_dseg_pf()
{return dseg_pf;}

double HXM::get_mflwseg_pf()
{return mflwseg_pf;}

double HXM::get_mflwvapseg_pf()
{return mflwvapseg_pf;}

double HXM::get_lcseg_fa()
{return 4.0*lseg_fa*thick_fa*memsur_corr/(2.0*lseg_fa*memsur_corr+2.0*thick_fa);}

double HXM::get_lcseg_pf()
{return 4.0*lseg_pf*thick_pf*memsur_corr/(2.0*lseg_pf*memsur_corr+2.0*thick_pf);}
