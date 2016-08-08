/* 
This is used to generate input and out files, written by Zhiming Gao (gaoz@ornl.gov)
*/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <new>
#include <typeinfo>
#include <exception>
#include <stdexcept>
#include <string>
#include <ios>
#include "iofile.h"

using namespace std;

namespace DeviceIO
{
	double length_dev, deep_dev, height_feedchanel, height_permeatechanel, thick_membrane;

	double poreradius_membrane, porosity_membrane, tortuosity_membrane, thermocond_membrane;

	double mass_feedflow, mass_permeateflow, temp_feedflow, rh_feedflow, press_feedflow, temp_permeateflow, rh_permeateflow, press_permeateflow; 

	int segment_length, segment_deep, layernum_membrane, flowtype_device;

	int flowcondition_permeate, membrane_model;

	double htc_correction, ftc_correction, memsurf_correction, mt_impact_ht;
	
	ofstream *output;
	
	string outputFilename;
}	
	
void DeviceIO::inputfile()
{
	ifstream *input;
	string ibuffer = "";
	string inputFilename = "";
	cout << endl << "Enter name of input file: " << endl;
	for (; ;)
	{
		cin >> inputFilename;
		input = new ifstream();
		input->open(inputFilename.c_str(), ios::in);
		if (input->is_open() == false)
		{
			cout << "File does not exist in the current directory." << endl;
			cout << "Enter new filename." << endl;
			delete input;
			input = 0;
		}
		else
		{break;}
	}
	//ofstream *output;
	string outputFilename = "";
	cout << endl << "Enter name of output file: " << endl;
	cin >> outputFilename;
	output = new ofstream();
	output->open(outputFilename.c_str(), ios::out);
	
	cout << endl << "Reading input data..." << endl;
	while (*input >> ibuffer)
	{
		if (ibuffer.compare("length_dev") == 0){*input >> length_dev;	continue;}
		if (ibuffer.compare("deep_dev") == 0){*input >> deep_dev;	continue;}
		if (ibuffer.compare("height_feedchanel") == 0){*input >> height_feedchanel;	continue;}
		if (ibuffer.compare("height_permeatechanel") == 0){*input >> height_permeatechanel;	continue;}
		if (ibuffer.compare("thick_membrane") == 0){*input >> thick_membrane;	continue;}
		if (ibuffer.compare("segment_length") == 0){*input >> segment_length;	continue;}
		if (ibuffer.compare("segment_deep") == 0){*input >> segment_deep;	continue;}
		if (ibuffer.compare("layernum_membrane") == 0){*input >> layernum_membrane;	continue;}
		if (ibuffer.compare("flowtype_device") == 0){*input >> flowtype_device;	continue;}
		if (ibuffer.compare("poreradius_membrane") == 0){*input >> poreradius_membrane;	continue;}
		if (ibuffer.compare("porosity_membrane") == 0){*input >> porosity_membrane;	continue;}
		if (ibuffer.compare("tortuosity_membrane") == 0){*input >> tortuosity_membrane;	continue;}
		if (ibuffer.compare("thermocond_membrane") == 0){*input >> thermocond_membrane;	continue;}
		if (ibuffer.compare("mass_feedflow") == 0){*input >> mass_feedflow;	continue;}
		if (ibuffer.compare("temp_feedflow") == 0){*input >> temp_feedflow;	continue;}
		if (ibuffer.compare("rh_feedflow") == 0){*input >> rh_feedflow;	continue;}
		if (ibuffer.compare("press_feedflow") == 0){*input >> press_feedflow;	continue;}
		if (ibuffer.compare("mass_permeateflow") == 0){*input >> mass_permeateflow;	continue;}
		if (ibuffer.compare("temp_permeateflow") == 0){*input >> temp_permeateflow;	continue;}
		if (ibuffer.compare("rh_permeateflow") == 0){*input >> rh_permeateflow;	continue;}
		if (ibuffer.compare("press_permeateflow") == 0){*input >> press_permeateflow;	continue;}
		if (ibuffer.compare("flowcondition_permeate") == 0){*input >> flowcondition_permeate;	continue;}
		if (ibuffer.compare("membrane_model") == 0){*input >> membrane_model;	continue;}
		if (ibuffer.compare("htc_correction") == 0){*input >> htc_correction;	continue;}
		if (ibuffer.compare("ftc_correction") == 0){*input >> ftc_correction;	continue;}
		if (ibuffer.compare("mt_impact_ht") == 0){*input >> mt_impact_ht;	continue;}
		if (ibuffer.compare("memsurf_correction") == 0){*input >> memsurf_correction;	continue;}
	}
	input->close();
	delete input;
	input = 0;
}


void DeviceIO::outputfile(int lnvalue, int dnvalue, int laynvalue, int flowtype, double tfa_in, double tfa_out, double humidfa_in, double humidfa_out, double rhfa_in, double rhfa_out, double Pvapfa_in, double Pvapfa_out,
	double tpf_in, double tpf_out, double humidpf_in, double humidpf_out, double rhpf_in, double rhpf_out, double Pvappf_in, double Pvappf_out,
	double * *map_temp_fa, double * *map_temp_fam, double * *map_vap_fa, double * *map_humid_fa, double * *map_flow_fa, double * *map_pvap_fa,
	double * *map_temp_pf, double * *map_temp_pfm, double * *map_vap_pf, double * *map_humid_pf, double * *map_flow_pf, double * *map_pvap_pf, double * *map_permeat_vap)
{
	if (flowtype == 1)
	{
		*output << "Notice: it is a parallel uncounterflow" << endl;
		*output << "tfa_in(c), tfa_out(c), humidfa_in(kg/kg), humidfa_out(kg/kg), rhfa_in, rhfa_out, Pvapfa_in(pa), Pvapfa_out(pa),	tpf_in(c), tpf_out(c), humidpf_in(kg/kg), humidpf_out(kg/kg), rhpf_in, rhpf_out, Pvappf_in(pa), Pvappf_out(pa)" << endl;
		*output << tfa_in - 273.15 << ", " << tfa_out - 273.15 << ", " << humidfa_in << ", " << humidfa_out << ", " << rhfa_in << ", " << rhfa_out << ", " << Pvapfa_in << ", " << Pvapfa_out << ", "
				<< tpf_in - 273.15 << ", " << tpf_out - 273.15 << ", " << humidpf_in << ", " << humidpf_out << ", " << rhpf_in << ", " << rhpf_out << ", " << Pvappf_in << ", " << Pvappf_out << endl;

		*output << "X/L, tfa(c), tfam(c), humidfa(kg/kg), pvapfa(pa), flfa(kg/s)@per_segment_per_channel, tpf(c), tpfm(c), humidpf(kg/kg), pvappf(pa), flpf(kg/s)@per_segment_per_channel, vappermeat(kg/m2/s)_from_feed-side_to_permeate-side@per_segment_per_channel" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i << ", " << map_temp_fa[0][i] - 273.15 << ", " << map_temp_fam[0][i] - 273.15 << ", " << map_humid_fa[0][i] << ", " << map_pvap_fa[0][i] << ", " << map_flow_fa[0][i] << ", "
				<< map_temp_pf[0][i] - 273.15 << ", " << map_temp_pfm[0][i] - 273.15 << ", " << map_humid_pf[0][i] << ", " << map_pvap_pf[0][i] << ", " << map_flow_pf[0][i] << ", " << map_permeat_vap[0][i] << endl;
		}  
	}
	else if (flowtype == 2)
	{
		*output << "Notice: it is a parallel counterflow" << endl;
		*output << "tfa_in(c), tfa_out(c), humidfa_in(kg/kg), humidfa_out(kg/kg), rhfa_in, rhfa_out, Pvapfa_in(pa), Pvapfa_out(pa),	tpf_in(c), tpf_out(c), humidpf_in(kg / kg), humidpf_out(kg / kg), rhpf_in, rhpf_out, Pvappf_in(pa), Pvappf_out(pa)" << endl;
		*output << tfa_in - 273.15 << ", " << tfa_out - 273.15 << ", " << humidfa_in << ", " << humidfa_out << ", " << rhfa_in << ", " << rhfa_out << ", " << Pvapfa_in << ", " << Pvapfa_out << ", "
				<< tpf_in - 273.15 << ", " << tpf_out - 273.15 << ", " << humidpf_in << ", " << humidpf_out << ", " << rhpf_in << ", " << rhpf_out << ", " << Pvappf_in << ", " << Pvappf_out << endl;

		*output << "X/L, tfa(c), tfam(c), humidfa(kg/kg), pvapfa(pa), flfa(kg/s)@per_segment_per_channel, tpf(c), tpfm(c), humidpf(kg/kg), pvappf(pa), flpf(kg/s)@per_segment_per_channel, vappermeat(kg/m2/s)_from_feed-side_to_permeate-side@per_segment_per_channel" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i << ", " << map_temp_fa[0][i] - 273.15 << ", " << map_temp_fam[0][i] - 273.15 << ", " << map_humid_fa[0][i] << ", " << map_pvap_fa[0][i] << ", " << map_flow_fa[0][i] << ", "
				<< map_temp_pf[0][i] - 273.15 << ", " << map_temp_pfm[0][i] - 273.15 << ", " << map_humid_pf[0][i] << ", " << map_pvap_pf[0][i] << ", " << map_flow_pf[0][i] << ", " << map_permeat_vap[0][i] << endl;
		}
	}
	else
	{
		*output << "Notice: it is a crossflow" << endl;
		*output << "tfa_in(c), tfa_out(c), humidfa_in(kg/kg), humidfa_out(kg/kg), rhfa_in, rhfa_out, Pvapfa_in(pa), Pvapfa_out(pa),	tpf_in(c), tpf_out(c), humidpf_in(kg / kg), humidpf_out(kg / kg), rhpf_in, rhpf_out, Pvappf_in(pa), Pvappf_out(pa)" << endl;
		*output << tfa_in - 273.15 << ", " << tfa_out - 273.15 << ", " << humidfa_in << ", " << humidfa_out << ", " << rhfa_in << ", " << rhfa_out << ", " << Pvapfa_in << ", " << Pvapfa_out << ", "
				<< tpf_in - 273.15 << ", " << tpf_out - 273.15 << ", " << humidpf_in << ", " << humidpf_out << ", " << rhpf_in << ", " << rhpf_out << ", " << Pvappf_in << ", " << Pvappf_out << endl;

		*output << "tfa(c)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_temp_fa[j][i] - 273.15); }
			*output << "" << endl;
		}

		*output << "tfam(c)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{
			*output << ", " << j;
		}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_temp_fam[j][i] - 273.15); }
			*output << "" << endl;
		}

		*output << "tpfm(c)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{
			*output << ", " << j;
		}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_temp_pfm[j][i] - 273.15); }
			*output << "" << endl;
		}

		*output << "tpf(c)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_temp_pf[j][i] - 273.15); }
			*output << "" << endl;
		}

		*output << "humidfa(kg/kg)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_humid_fa[j][i]); }
			*output << "" << endl;
		}

		*output << "humidpf(kg/kg)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << (map_humid_pf[j][i]); }
			*output << "" << endl;
		}

		*output << "pvapfa(pa)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{
			*output << ", " << j;
		}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << map_pvap_fa[j][i]; }
			*output << "" << endl;
		}

		*output << "pvappf(kg/kg)_MAP" << endl;
		*output << "D/L_index";
		for (int j = 0; j < lnvalue; j++)
		{
			*output << ", " << j;
		}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << map_pvap_pf[j][i]; }
			*output << "" << endl;
		}

		*output << "flfa(kg/s)_MAP @ per segment per channel" << endl;
		*output << "D/L index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << map_flow_fa[j][i]; }
			*output << "" << endl;
		}

		*output << "flpf(kg/s)_MAP @ per segment per channel" << endl; 
		*output << "D/L index";
		for (int j = 0; j < lnvalue; j++)
		{*output << ", " << j;}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << map_flow_pf[j][i]; }
			*output << "" << endl;
		}

		*output << "vappermeat(kg/m2/s)_MAP from feed-side to permeate-side @ per segment per channel" << endl; 
		*output << "D/L index";
		for (int j = 0; j < lnvalue; j++)
		{
			*output << ", " << j;
		}
		*output << "" << endl;
		for (int i = 0; i < dnvalue; i++)
		{
			*output << i;
			for (int j = 0; j < lnvalue; j++) { *output << ", " << map_permeat_vap[j][i]; }
			*output << "" << endl;
		}
	}

	output->close();
	delete output;
	output = 0;
}
