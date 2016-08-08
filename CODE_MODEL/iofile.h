#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#include <string>
using namespace std;

namespace DeviceIO
{
	void inputfile();
	// Gets data from a text file

	void outputfile(int lnvalue, int dnvalue, int laynvalue, int ftype, double tfa_in, double tfa_out, double humidfa_in, double humidfa_out, double rhfa_in, double rhfa_out, double Pvapfa_in, double Pvapfa_out,
		double tpf_in, double tpf_out, double humidpf_in, double humidpf_out, double rhpf_in, double rhpf_out, double Pvappf_in, double Pvappf_out,
		double * *map_temp_fa, double * *map_temp_fam, double * *map_vap_fa, double * *map_humid_fa, double * *map_flow_fa, double * *map_pvap_fa,
		double * *map_temp_pf, double * *map_temp_pfm, double * *map_vap_pf, double * *map_humid_pf, double * *map_flow_pf, double * *map_pvap_pf, double * *map_permeat_vap);
	// output data to a text file

	extern double length_dev, deep_dev, height_feedchanel, height_permeatechanel, thick_membrane;

	extern double poreradius_membrane, porosity_membrane, tortuosity_membrane, thermocond_membrane;

	extern double mass_feedflow, mass_permeateflow, temp_feedflow, rh_feedflow, press_feedflow, temp_permeateflow, rh_permeateflow, press_permeateflow;

	extern int segment_length, segment_deep, layernum_membrane, flowtype_device;

	extern int flowcondition_permeate, membrane_model;

	extern double htc_correction, ftc_correction, memsurf_correction, mt_impact_ht;

	extern ofstream *output;
}
#endif