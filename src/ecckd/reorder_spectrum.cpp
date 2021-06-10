// reorder_spectrum.cpp - Reorder spectrum in each band in order of increasing absorption
//
// Copyright (C) 2020- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.
//
// Author:  Robin Hogan
// Email:   r.j.hogan@ecmwf.int

#include <vector>
#include <string>
#include <algorithm>

#include "DataFile.h"
#include "radiative_transfer_lw.h"
#include "radiative_transfer_sw.h"
#include "read_spectrum.h"
#include "planck_function.h"
#include "heating_rate.h"
#include "write_order.h"
#include "calc_cost_function_sw.h"

/// Comparison structure for rearranging a std::vector of indices so
/// that they index an adept::Vector of data in ascending order
struct MyCompare {
  MyCompare(const adept::Vector& data) : data_(data) {}
  bool operator()(int i1, int i2) { return data_[i1] < data_[i2]; }
  const adept::Vector& data_;
};

int
main(int argc, const char* argv[])
{
  using namespace adept;

  // Names of input and output files
  std::string input, output, ssi_file_name;

  // Threshold optical depth
  Real threshold_optical_depth = 0.5;

  // CONFIGURATION

  // Read configuration information from command-line and first file
  // on command-line
  DataFile config(argc, argv);

  std::string log_level;
  if (config.read(log_level, "log_level")) {
    set_log_level(log_level);
  }

  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  bool do_sw = false;
  if (config.read(ssi_file_name, "ssi")) {
    do_sw = true;
    LOG << "Assuming shortwave spectral region (ssi provided)\n";
  }
  else {
    LOG << "Assuming longwave spectral region (ssi not provided)\n";
  }

  std::string molecule;

  // Index of profile to extract from file
  int iprofile = 0;
  config.read(iprofile, "iprofile");

  config.read(threshold_optical_depth, "threshold_optical_depth");

  // READ SPECTRAL OPTICAL DEPTH

  LOG << "Reading " << input << "\n";

  // Wavenumber and wavenumber spacing, cm-1
  Vector wavenumber_cm_1, d_wavenumber_cm_1;
  // Half-level pressure, Pa
  Vector pressure_hl;
  Matrix optical_depth;
  Vector temp_hl; // Temperature, not used
  Real reference_surface_vmr;
  Vector vmr_fl; // Not used

  read_spectrum(input, iprofile, pressure_hl, temp_hl, wavenumber_cm_1, 
		d_wavenumber_cm_1, optical_depth, molecule,
		reference_surface_vmr, vmr_fl);

  // Override molecule name if present
  config.read(molecule, "molecule");

  int nlay = pressure_hl.size()-1;
  int nwav = wavenumber_cm_1.size();

  LOG << nlay << " layers\n";
  LOG << nwav << " spectral points\n";

  Matrix flux_dn(nlay+1,nwav), flux_up(nlay+1,nwav);

  if (!do_sw) {

    // COMPUTE IDEALIZED PLANCK FUNCTION
  
    LOG << "Computing Planck function\n";
  
    // Idealized linearly increasing temperature profile with
    // log(pressure)
    Vector log_p_interp = {log(1.0), log(100000.0)};
    Vector temp_interp  = {273.15-100.0, 273.15+15.0};
  
    Vector temperature_hl = interp(log_p_interp, temp_interp, log(pressure_hl));
    Matrix planck_hl(nlay+1,nwav);

    planck_function(temperature_hl, wavenumber_cm_1, d_wavenumber_cm_1,
		    planck_hl);

    Vector surf_planck(nwav);
    planck_function(temperature_hl(end), wavenumber_cm_1, d_wavenumber_cm_1, 
		    surf_planck);

    // RADIATIVE TRANSFER

    LOG << "Performing longwave radiative transfer\n";
    
    Vector surf_emissivity(nwav);
    surf_emissivity = 1.0;
    
    radiative_transfer_lw(planck_hl, optical_depth, surf_emissivity, surf_planck,
			  flux_dn, flux_up);
    planck_hl.clear();

  }
  else {

    Real cos_sza = REFERENCE_COS_SZA;

    Vector ssi;

    LOG << "Reading " << ssi_file_name << "\n";
    DataFile ssi_file(ssi_file_name);
    ssi_file.read(ssi, "solar_spectral_irradiance");

    radiative_transfer_direct_sw(cos_sza, ssi,
				 optical_depth, flux_dn);
    flux_up = 0.0;

  }

  Vector column_optical_depth = sum(optical_depth, 0);
  //  optical_depth.clear();
  
  // CALCULATE HEATING RATE

  LOG << "Computing heating rate\n";

  Matrix hr(nlay,nwav);
  heating_rate(pressure_hl, flux_dn, flux_up, hr);

  if (!do_sw) {
    // We are only interested in cooling
    LOG << "Locating peak cooling\n";
    hr.where(hr > 0.0) = 0.0;
  }

  Vector pseudo_height = log(pressure_hl(end)) - 0.5*( log(pressure_hl(range(0,end-1)))
						      +log(pressure_hl(range(1,end))));
  Vector d_height = log(pressure_hl(range(1,end))) - log(pressure_hl(range(0,end-1)));
  // This is actually the peak heating in the shortwave
  Vector peak_cooling_height = sum(hr * spread<1>(d_height*pseudo_height,nwav), 0)
    / sum(hr * spread<1>(d_height,nwav), 0);

  // If column optical depth is less than a threshold, don't sort by
  // the heating rate pressure, but by the optical depth itself
  if (threshold_optical_depth > 0.0) {
    peak_cooling_height.where(column_optical_depth < threshold_optical_depth)
      = -threshold_optical_depth + column_optical_depth;
  }

  // Find the pseudo height at which optical depth down from TOA
  // reaches threshold
  LOG << "Finding height at which optical depth from TOA reaches "
	<< threshold_optical_depth << "\n";

  Vector pseudo_height_hl = log(pressure_hl(end)) - log(pressure_hl);
  Vector od_threshold_height(nwav);
  for (int iwav = 0; iwav < nwav; ++iwav) {
    if (column_optical_depth(iwav) <= threshold_optical_depth) {
      od_threshold_height(iwav) = column_optical_depth(iwav) - threshold_optical_depth;
    }
    else {
      Real cum_od = 0.0;
      od_threshold_height(iwav) = 0.0;
      for (int ilay = 0; ilay < nlay; ++ilay) {
	Real next_cum_od = cum_od + optical_depth(ilay,iwav);
	if (next_cum_od >= threshold_optical_depth) {
	  // Enclosed point where optical depth exceeds threshold
	  od_threshold_height(iwav) 
	    = ((threshold_optical_depth-cum_od)*pseudo_height_hl(ilay+1)
	       + (next_cum_od-threshold_optical_depth)*pseudo_height_hl(ilay))
	    / std::max(1.0e-12,optical_depth(ilay,iwav));
	  if (od_threshold_height(iwav) > 30.0) {
	    throw;
	  }
	  break;
	}
	cum_od = next_cum_od;
      }
    }
  }

  if (do_sw) {
    LOG << "Using height at which optical depth from TOA reaches "
	<< threshold_optical_depth << "\n";
    peak_cooling_height = od_threshold_height;
  }

  Vector band_bound1, band_bound2;
  int nband = 1;
  if (config.exist("wavenumber1")) {
    config.read(band_bound1, "wavenumber1");
    config.read(band_bound2, "wavenumber2");
    nband = band_bound1.size();
  }
  else {
    band_bound1.resize(1);
    band_bound2.resize(1);
    band_bound1(0) = std::max(0.0, wavenumber_cm_1(0)-d_wavenumber_cm_1(0));
    band_bound2(0) = wavenumber_cm_1(end)+d_wavenumber_cm_1(end);
  }

  if (nband <= 0) {
    ERROR << "Failure to interpret wavenumber1 and wavenumber2 as a list of band boundaries";
    THROW(PARAMETER_ERROR);
  }
  if (nband == 1) {
    LOG << "Treating the entire spectrum as one band\n";
  }
  else {
    LOG << "Splitting the spectrum into " << nband << " bands\n";
  }

  if (!do_sw) {
    LOG << "Sorting by peak cooling\n";
  }
  else {
    LOG << "Sorting by peak heating\n";
  }

  std::vector<int> g_index(nwav);
  for (int jwav = 0; jwav < nwav; ++jwav) {
    g_index[jwav] = jwav;
  }
  MyCompare my_compare(peak_cooling_height);

  // Bounds that clamp to the range of the data
  Vector band_bound_clamp1, band_bound_clamp2;
  band_bound_clamp1 = band_bound1;
  band_bound_clamp2 = band_bound2;
  band_bound_clamp1(0)   = std::max(wavenumber_cm_1(0),band_bound1(0));
  band_bound_clamp2(end) = std::min(wavenumber_cm_1(end),band_bound2(end));

  intVector iband(nwav);
  iband = -1;
  for (int jband = 0; jband < nband; ++jband) {
    LOG << "  Band " << jband << ": " << band_bound_clamp1(jband)
	<< "-" << band_bound_clamp2(jband) << " cm-1\n";
    intVector index;
    if (jband < nband-1) {
      index = find(   wavenumber_cm_1 >= band_bound1(jband)
		   && wavenumber_cm_1 <  band_bound2(jband));
    }
    else {
      index = find(   wavenumber_cm_1 >= band_bound1(jband)
		   && wavenumber_cm_1 <= band_bound2(jband));
    }
    iband(index) = jband;
    int index1 = index(0);
    int index2 = index(end);
    std::vector<int>::iterator start  = g_index.begin() + index1;
    std::vector<int>::iterator finish = g_index.begin() + index2 + 1;
    std::stable_sort(start, finish, my_compare);
  }

  intVector ordered_index(&g_index[0], dimensions(nwav)); // Point to data in g_index

  intVector rank(nwav);
  rank(ordered_index) = range(0,nwav-1);

  LOG << "Writing " << output << "\n";

  // Get configuration information as a string
  std::string config_str;
  config.read(config_str);  

  write_order(output, argc, argv, molecule, config_str,
	      band_bound_clamp1, band_bound_clamp2,
	      wavenumber_cm_1, d_wavenumber_cm_1,
	      iband, rank, ordered_index, column_optical_depth, peak_cooling_height);

}
