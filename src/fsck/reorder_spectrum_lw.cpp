#include <vector>
#include <string>
#include <algorithm>

#include "DataFile.h"
#include "radiative_transfer_lw.h"
#include "read_spectrum.h"
#include "planck_function.h"
#include "heating_rate.h"
#include "write_order.h"

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
  std::string input, output;

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
  Matrix flux_dn(nlay+1,nwav), flux_up(nlay+1,nwav);

  radiative_transfer_lw(planck_hl, optical_depth, surf_emissivity, surf_planck,
			flux_dn, flux_up);
  Vector column_optical_depth = sum(optical_depth, 0);
  optical_depth.clear();
  planck_hl.clear();
  
  // CALCULATE HEATING RATE

  LOG << "Computing heating rate\n";

  Matrix hr(nlay,nwav);
  heating_rate(pressure_hl, flux_dn, flux_up, hr);

  // We are only interested in cooling
  LOG << "Locating peak cooling\n";
  hr.where(hr > 0.0) = 0.0;

  Vector pseudo_height = log(pressure_hl(end)) - 0.5*( log(pressure_hl(range(0,end-1)))
						      +log(pressure_hl(range(1,end))));
  Vector d_height = log(pressure_hl(range(1,end))) - log(pressure_hl(range(0,end-1)));
  Vector peak_cooling_height = sum(hr * spread<1>(d_height*pseudo_height,nwav), 0)
    / sum(hr * spread<1>(d_height,nwav), 0);

  // If column optical depth is less than a threshold, don't sort by
  // the heating rate pressure, but by the optical depth itself
  if (threshold_optical_depth > 0.0) {
    peak_cooling_height.where(column_optical_depth < threshold_optical_depth)
      = -threshold_optical_depth + column_optical_depth;
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

  LOG << "Sorting by peak cooling\n";

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
    intVector index = find(   wavenumber_cm_1 >= band_bound1(jband)
			   && wavenumber_cm_1 <  band_bound2(jband));
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
