#include <string>

#include "read_merged_spectrum.h"
#include "planck_function.h"
#include "radiative_transfer_lw.h"
#include "radiative_transfer_sw.h"
#include "heating_rate.h"
#include "calc_cost_function_lw.h"
#include "calc_cost_function_sw.h"
#include "OutputDataFile.h"
#include "single_gas_data.h"
#include "equipartition.h"
#include "write_standard_attributes.h"
#include "cumsum.h"

using namespace adept;


// Return the Planck-weighted-median of the sorting variable indexed
// by "index", where "rank" provides the ordering of sorting_variable
Real
calc_median_sorting_variable(const Vector& sorting_variable,
			     const Vector& planck,
			     int i1, int i2) {
  Real half_planck = 0.5 * sum(planck(range(i1,i2)));
  Real cum_planck = 0.0;
  int iind = i1;
  for ( ; iind < i2; ++iind) {
    cum_planck += planck(iind);
    if (cum_planck >= half_planck) {
      break;
    }
  }
  return sorting_variable(iind);
}

// Calculate the "fitted" optical depth averaged over a range of
// wavenumbers from indices i1 to i2, using the Planck function to
// weight at each height
Vector fit_optical_depth_lw(const std::string& averaging_method,
			    int i1, int i2,
			    const Matrix& planck_hl,
			    const Matrix& metric) {
  Vector optical_depth_fit;

  if (averaging_method == "linear") {
    optical_depth_fit = sum(metric(__,range(i1,i2)) * planck_hl(range(1,end),range(i1,i2)), 1)
      / sum(planck_hl(range(1,end),range(i1,i2)),1);
  }
  else if (averaging_method == "transmission") {
    optical_depth_fit = min(0.9999999999999999,sum(metric(__,range(i1,i2)) * planck_hl(range(1,end),range(i1,i2)), 1)
			    / sum(planck_hl(range(1,end),range(i1,i2)),1));
    optical_depth_fit = abs(-log(1.0-optical_depth_fit)/LW_DIFFUSIVITY);
  }
  else if (averaging_method == "transmission-2") {
    optical_depth_fit = min(0.9999999999999999,sum(metric(__,range(i1,i2)) * planck_hl(range(1,end),range(i1,i2)), 1)
			    / sum(planck_hl(range(1,end),range(i1,i2)),1));
    optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*2.0));
  }
  else if (averaging_method == "square-root") {
    optical_depth_fit = sum(metric(__,range(i1,i2)) * planck_hl(range(1,end),range(i1,i2)), 1)
      / sum(planck_hl(range(1,end),range(i1,i2)),1);
    optical_depth_fit *= optical_depth_fit;
  }
  else if (averaging_method == "logarithmic") {
    optical_depth_fit.resize(metric.size(0));
    for (int iz = 0; iz < metric.size(0); ++iz) {
      intVector iindex = i1+find(metric(iz,range(i1,i2)) > 0.0);
      Vector od_nonzero = metric(iz,iindex);
      if (od_nonzero.size() == i2-i1+1) {
	// Pure logarithmic average
	optical_depth_fit(iz) = exp(sum(log(od_nonzero)*planck_hl(iz+1,iindex))
				    / sum(planck_hl(iz,iindex)));
      }
      else if (od_nonzero.empty()) {
	// No non-zero data
	optical_depth_fit(iz) = 0.0;
      }
      else {
	// Some zeros: logarithmic average of non-zeros then linear average with zeros
	optical_depth_fit(iz) = exp(sum(log(od_nonzero)*planck_hl(iz+1,iindex))
				    / sum(planck_hl(iz,iindex)))
	  * (static_cast<Real>(iindex.size())/static_cast<Real>(i2-i1+1));
      }
    }
  }
  else {
    ERROR << "Averaging method \"" << averaging_method << "\" not understood";
    THROW(PARAMETER_ERROR);
  }
  return optical_depth_fit;
}


// Calculate the "fitted" optical depth averaged over a range of
// wavenumbers from indices i1 to i2, using the spectral solar
// irradiance
Vector fit_optical_depth_sw(const std::string& averaging_method,
			    int i1, int i2,
			    const Vector& ssi,
			    const Matrix& metric) {
  Vector optical_depth_fit;
  int nz = metric.size(0);
  Real norm_factor = 1.0 / sum(ssi(range(i1,i2)));
  if (averaging_method == "linear") {
    optical_depth_fit = sum(metric(__,range(i1,i2)) * spread<0>(ssi(range(i1,i2)),nz),1) * norm_factor;
  }
  else if (averaging_method == "transmission") {
    optical_depth_fit = min(0.9999999999999999,sum(metric(__,range(i1,i2)) * spread<0>(ssi(range(i1,i2)),nz),1))
      * norm_factor;
    optical_depth_fit = abs(-log(1.0-optical_depth_fit)/LW_DIFFUSIVITY);
  }
  else if (averaging_method == "transmission-2") {
    optical_depth_fit = min(0.9999999999999999,sum(metric(__,range(i1,i2)) * spread<0>(ssi(range(i1,i2)),nz),1))
      * norm_factor;
    optical_depth_fit = abs(-log(1.0-optical_depth_fit)/(LW_DIFFUSIVITY*2.0));
  }
  else if (averaging_method == "square-root") {
    optical_depth_fit = sum(metric(__,range(i1,i2)) * spread<0>(ssi(range(i1,i2)),nz),1)
      * norm_factor;
    optical_depth_fit *= optical_depth_fit;
  }
  else if (averaging_method == "logarithmic"
	   || averaging_method == "total-transmission") {
    optical_depth_fit.resize(nz);
    for (int iz = 0; iz < nz; ++iz) {
      intVector iindex = i1+find(metric(iz,range(i1,i2)) > 0.0);
      Vector od_nonzero = metric(iz,iindex);
      if (od_nonzero.size() == i2-i1+1) {
	// Pure logarithmic average
	optical_depth_fit(iz) = exp(sum(log(od_nonzero)*ssi(iindex))
				    / sum(ssi(iindex)));
      }
      else if (od_nonzero.empty()) {
	// No non-zero data
	optical_depth_fit(iz) = 0.0;
      }
      else {
	// Some zeros: logarithmic average of non-zeros then linear average with zeros
	optical_depth_fit(iz) = exp(sum(log(od_nonzero)*ssi(iindex))
				    / sum(ssi(iindex)))
	  * (static_cast<Real>(iindex.size())/static_cast<Real>(i2-i1+1));
      }
    }
  }
  else {
    ERROR << "Averaging method \"" << averaging_method << "\" not understood";
    THROW(PARAMETER_ERROR);
  }
  return optical_depth_fit;
}


// Calculate the "fitted" optical depth averaged over a range of
// wavenumbers from indices i1 to i2, using the spectral solar
// irradiance and the total-transmission shortwave method
Vector fit_optical_depth_sw_total_trans(
			    int i1, int i2,
			    const Vector& ssi,
			    const Matrix& bg_od, const Matrix& od) {
  int nz = od.size(0);
  Vector optical_depth_fit(nz);
  Vector flux_dn; flux_dn = ssi(range(i1,i2));
  Vector bg_flux_dn; bg_flux_dn = flux_dn;
  Real bb_flux_dn_top, bb_flux_dn_base;
  Real bb_bg_flux_dn_top, bb_bg_flux_dn_base;

  bb_flux_dn_top = bb_bg_flux_dn_top = sum(flux_dn);

  Real norm_factor = 1.0 / sum(ssi(range(i1,i2)));
  Real bg_od_fit;

  for (int iz = 0; iz < nz; ++iz) {
    // Factor of 2.0 is because we consider a 60 degree solar zenith angle
    bg_flux_dn *= exp(-2.0 * bg_od(iz,range(i1,i2)));
    flux_dn *= exp(-2.0 * (bg_od(iz,range(i1,i2))+od(iz,range(i1,i2))));
    bb_bg_flux_dn_base = sum(bg_flux_dn);
    bb_flux_dn_base = sum(flux_dn);
    if (bb_bg_flux_dn_base > 0.0 && bb_flux_dn_base > 0.0) {
      bg_od_fit = -0.5*log(bb_bg_flux_dn_base / bb_bg_flux_dn_top);
      optical_depth_fit(iz) = -0.5*log(bb_flux_dn_base / bb_flux_dn_top) - bg_od_fit;
    }
    else {
      optical_depth_fit = sum(od(__,range(i1,i2)) * spread<0>(ssi(range(i1,i2)),nz),1) * norm_factor;
    }
    bb_flux_dn_top = bb_flux_dn_base;
    bb_bg_flux_dn_top = bb_bg_flux_dn_base;
  }
  return optical_depth_fit;
}

class CkdEquipartition : public Equipartition {
public:
  CkdEquipartition() { }
  void init_lw(std::string am, Real fw, const Vector& lw,
	       const Vector& prhl, const Vector& se,
	       const Vector& sp, const Vector& fds, const Vector& fut,
	       const Matrix& plhl, const Matrix& bod, const Matrix& met,
	       const Matrix& h, int i1, int i2) {
    do_sw = false;
    averaging_method = am;
    flux_weight = fw;
    layer_weight = lw;
    pressure_hl = prhl;
    surf_emissivity = se(range(i1,i2));
    surf_planck = sp(range(i1,i2));
    flux_dn_surf = fds(range(i1,i2));
    flux_up_toa = fut(range(i1,i2));
    planck_hl = plhl(__,range(i1,i2));
    bg_optical_depth = bod(__,range(i1,i2));
    metric = met(__,range(i1,i2));
    hr = h(__,range(i1,i2));
    npoints = i2-i1+1;
    total_comp_cost = 0.0;
    
    set_resolution(1.0 / npoints);
    set_parallel(true);
    set_verbose(true);
    set_minimize_frac_range(true);
    debug_partition = false;
  }

  void init_sw(std::string am, Real fw, const Vector& lw, Real cs,
	       const Vector& prhl, const Vector& si, Real sa,
	       const Vector& fds, const Vector& fut, 
	       const Matrix& bod, const Matrix& met,
	       const Matrix& h, int i1, int i2) {
    do_sw = true;
    averaging_method = am;
    flux_weight = fw;
    layer_weight = lw;
    cos_sza = cs;
    pressure_hl = prhl;
    ssi = si(range(i1,i2));
    surf_albedo = sa;
    flux_dn_surf = fds(range(i1,i2));
    flux_up_toa = fut(range(i1,i2));
    bg_optical_depth = bod(__,range(i1,i2));
    metric = met(__,range(i1,i2));
    hr = h(__,range(i1,i2));
    npoints = i2-i1+1;
    total_comp_cost = 0.0;
    
    set_resolution(1.0 / npoints);
    set_parallel(true);
    set_verbose(true);
    set_minimize_frac_range(true);
    debug_partition = false;

  }

  void init_sw_extras(const Vector& fdsl, const Vector& futl,
		      const Vector& fdsh, const Vector& futh,
		      Real mins, Real maxs,
		      const Matrix& hl, const Matrix& hh,
		      int i1, int i2) {
    flux_dn_surf_low  = fdsl(range(i1,i2));
    flux_up_toa_low   = futl(range(i1,i2));
    flux_dn_surf_high = fdsh(range(i1,i2));
    flux_up_toa_high  = futh(range(i1,i2));
    min_scaling = mins;
    max_scaling = maxs;
    hr_low  = hl(__,range(i1,i2));
    hr_high = hh(__,range(i1,i2));
  }


  int lower_index(ep_real bound) const {
    return static_cast<int>(std::ceil(bound*(npoints-1)));
  }
  int upper_index(ep_real bound) const {
    return static_cast<int>(std::floor(bound*(npoints-1)));
  }

  // Need to take care that this is thread safe: soft_link() ensures
  // no reference counting for read-only objects
  ep_real calc_error(ep_real bound1, ep_real bound2) {
    std::cout << "." << std::flush;
    int ibound1 = lower_index(bound1);
    //int ibound2 = std::max(ibound1,upper_index(bound2));
    int ibound2 = upper_index(bound2);
    //std::cout << "[" << bound1 << " " << bound2 << " " << ibound1 << " " << ibound2 << "]" << std::flush;
    if (ibound1 < 0 || ibound2 >= npoints) {
      std::cout << "*** ERROR: requested bounds " << bound1 << "-" << bound2
		<< " corresponding to indices " << ibound1 << "-" << ibound2 
		<< " outside valid range 0-" << npoints-1 << std::endl;
      throw(PROCESSING_ERROR);
    }
    else if (ibound2 < ibound1) {
      std::cout << "*** ERROR: requested indices out of order: " 
		<< ibound1 << "-" << ibound2 << std::endl;
      throw(PROCESSING_ERROR);
    }

    //total_comp_cost += ibound2-ibound1+1;
    total_comp_cost += bound2-bound1;

    if (!do_sw) {
      Vector optical_depth_fit = fit_optical_depth_lw(averaging_method,
						      ibound1, ibound2,
						      planck_hl.soft_link(),
						      metric.soft_link());
	  
      return calc_cost_function_lw(pressure_hl.soft_link(),
				   planck_hl.soft_link()(__,range(ibound1,ibound2)),
				   surf_emissivity.soft_link()(range(ibound1,ibound2)),
				   surf_planck.soft_link()(range(ibound1,ibound2)),
				   bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
				   optical_depth_fit.soft_link(),
				   flux_dn_surf.soft_link()(range(ibound1,ibound2)),
				   flux_up_toa.soft_link()(range(ibound1,ibound2)),
				   hr.soft_link()(__,range(ibound1,ibound2)),
				   flux_weight, layer_weight.soft_link());
    }
    else {
      if (averaging_method == "total-transmission") {
	Vector optical_depth_fit = fit_optical_depth_sw_total_trans(
							ibound1, ibound2,
							ssi.soft_link(),
							bg_optical_depth.soft_link(),
							metric.soft_link());
	if (debug_partition) {
	  std::cerr << "  debug_partition_LOW\n";
	}
	Real cf_low = calc_cost_function_sw(cos_sza,
				     pressure_hl.soft_link(),
				     ssi.soft_link()(range(ibound1,ibound2)), surf_albedo,
				     bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
				     optical_depth_fit * min_scaling,
				     flux_dn_surf_low.soft_link()(range(ibound1,ibound2)),
				     flux_up_toa_low.soft_link()(range(ibound1,ibound2)),
				     hr_low.soft_link()(__,range(ibound1,ibound2)),
				     flux_weight, layer_weight.soft_link(),
				     intVector(), debug_partition);
	if (debug_partition) {
	  std::cerr << "  debug_partition_HIGH\n";
	}
	Real cf_high = calc_cost_function_sw(cos_sza,
				     pressure_hl.soft_link(),
				     ssi.soft_link()(range(ibound1,ibound2)), surf_albedo,
				     bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
				     optical_depth_fit * max_scaling,
				     flux_dn_surf_high.soft_link()(range(ibound1,ibound2)),
				     flux_up_toa_high.soft_link()(range(ibound1,ibound2)),
				     hr_high.soft_link()(__,range(ibound1,ibound2)),
				     flux_weight, layer_weight.soft_link(),
				     intVector(), debug_partition);
	if (debug_partition) {
	  std::cerr << "  debug_partition_MID\n";
	  Real cf = calc_cost_function_sw(cos_sza,
					  pressure_hl.soft_link(),
					  ssi.soft_link()(range(ibound1,ibound2)), surf_albedo,
					  bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
					  optical_depth_fit,
					  flux_dn_surf.soft_link()(range(ibound1,ibound2)),
					  flux_up_toa.soft_link()(range(ibound1,ibound2)),
					  hr.soft_link()(__,range(ibound1,ibound2)),
					  flux_weight, layer_weight.soft_link(),
					  intVector(), debug_partition);
	  //	  std::cout << "*** cost functions = " << cf_low << " " << cf << " " << cf_high << "\n";
	}
	return 0.5 * (cf_low + cf_high);
      }
      else {
	Vector optical_depth_fit = fit_optical_depth_sw(averaging_method,
							ibound1, ibound2,
							ssi.soft_link(),
							metric.soft_link());
	  
	return calc_cost_function_sw(cos_sza,
				     pressure_hl.soft_link(),
				     ssi.soft_link()(range(ibound1,ibound2)), surf_albedo,
				     bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
				     optical_depth_fit.soft_link(),
				     flux_dn_surf.soft_link()(range(ibound1,ibound2)),
				     flux_up_toa.soft_link()(range(ibound1,ibound2)),
				     hr.soft_link()(__,range(ibound1,ibound2)),
				     flux_weight, layer_weight.soft_link());
      }
    }
  }

  std::string averaging_method;
  Real flux_weight;
  Real cos_sza;
  Vector layer_weight;
  Vector pressure_hl;
  Vector ssi;
  Vector surf_emissivity, surf_planck;
  Real surf_albedo;
  Vector flux_dn_surf, flux_up_toa;
  Matrix planck_hl, bg_optical_depth, metric, hr;
  int npoints;
  ep_real total_comp_cost;
  bool do_sw;

  Vector flux_dn_surf_low, flux_up_toa_low;
  Vector flux_dn_surf_high, flux_up_toa_high;
  Real min_scaling, max_scaling;
  Matrix hr_low, hr_high;
  bool debug_partition = false;
};


// Main program 
int
main(int argc, const char* argv[])
{

  adept::set_array_print_style(PRINT_STYLE_MATLAB);

  // CONFIGURATION

  // Read configuration information from command-line and first file
  // on command-line
  DataFile config(argc, argv);

  std::string log_level;
  if (config.read(log_level, "log_level")) {
    set_log_level(log_level);
  }

  // Modify search path
  std::string mod_path;
  if (config.read(mod_path, "prepend_path")) {
    prepend_search_directory(mod_path);
  }
  if (config.read(mod_path, "append_path")) {
    append_search_directory(mod_path);
  }

  // Allow debugging
  //  set_trace_exceptions(true);

  // Name of output file
  std::string output, ssi_file_name;

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  bool do_sw = false;
  Real cos_sza = REFERENCE_COS_SZA;
  Real reference_albedo = 0.15;
  Vector ssi;

  bool debug_partition = false;
  config.read(debug_partition, "debug_partition");

  if (config.read(ssi_file_name, "ssi")) {
    do_sw = true;
    LOG << "Assuming shortwave spectral region (ssi provided)\n";
    LOG << "Reading " << ssi_file_name << "\n";
    DataFile ssi_file(ssi_file_name);
    ssi_file.read(ssi, "solar_spectral_irradiance");    
  }
  else {
    LOG << "Assuming longwave spectral region (ssi not provided)\n";
  }

  std::string molecule;

  // Index of profiles to extract from files
  int iprofile = 0;
  config.read(iprofile, "iprofile");

  Vector heating_rate_tolerance_in; // K d-1
  if (!config.read(heating_rate_tolerance_in, "heating_rate_tolerance")) {
    ERROR << "heating_rate_tolerance not defined";
    THROW(PARAMETER_ERROR);
  }

  // Aim for errors in each g point to be within 2% of each other
  Real tolerance_tolerance = 0.02;
  config.read(tolerance_tolerance, "tolerance_tolerance");
  //  Real heating_rate_upper_tolerance = heating_rate_tolerance * (1.0+tolerance_tolerance);
  //  Real heating_rate_lower_tolerance = heating_rate_tolerance * (1.0-tolerance_tolerance);

  int max_iterations = 60;
  config.read(max_iterations, "max_iterations");

  std::string averaging_method = "linear";
  config.read(averaging_method, "averaging_method");

  // Hogan (2012): 0.02 (K day-1 W-1 m2)^2
  Real flux_weight = 0.02;
  config.read(flux_weight, "flux_weight");

  Real repartition_factor = 1.0;
  config.read(repartition_factor, "repartition_factor");

  int repartition_repeat = 1;
  config.read(repartition_repeat, "repartition_repeat");

  // Wavenumber above which is it not safe to neglect Rayleigh
  // scattering in calculating upwelling fluxes
  Real max_no_rayleigh_wavenumber = 10000.0;
  config.read(max_no_rayleigh_wavenumber, "max_no_rayleigh_wavenumber");

  int ngas = 0;
  int nband = 0;
  std::string gas_str;

  Matrix planck_hl;
  Vector surf_planck;
  Vector band_bound1, band_bound2;

  // Surface albedo at every high-res wavenumber and just in the bands
  Vector albedo, band_albedo;

  std::vector<SingleGasData> single_gas_data;

  int nwav = 0; // Number of wavenumbers
  // Wavenumber, cm-1
  Vector wavenumber_cm_1;

  // If we consider clouds do that first
  std::string cloud_str;
  if (config.read(cloud_str, "cloud")) {
    LOG << "*** FINDING G POINTS FOR " << cloud_str << "\n";
    if (!do_sw) {
      ERROR << "Don't yet know how to sort cloud properties in the longwave";
      THROW(PARAMETER_ERROR);
    }

    // READ ORDERING
    std::string reordering_input;
    if (!config.read(reordering_input, cloud_str, "reordering_input")) {
      ERROR << "No reordering_input found";
      THROW(PARAMETER_ERROR);
    }
    LOG << "Reading " << reordering_input << "\n";
    DataFile order_file(reordering_input);
    Vector sorting_variable;
    intVector irank, iband;
    order_file.read(irank, "rank");
    order_file.read(iband, "band_number");
    order_file.read(band_bound1, "wavenumber1_band");
    order_file.read(band_bound2, "wavenumber2_band");
    order_file.read(sorting_variable, "sorting_variable");
    nband = band_bound1.size();

    // What is the maximum range of reflectance permitted for one g
    // point?
    Real max_reflectance_range = 0.26;
    config.read(max_reflectance_range, cloud_str, "max_reflectance_range");

    intVector n_g_points(nband);
    int ng = 0;

    std::vector<int> n_g_points_per_band;
    // Lower and upper g
    std::vector<int> rank1_per_g_point, rank2_per_g_point, band_num;
    // Heating rate error, K d-1
    std::vector<Real> error_per_g_point;
    // Median of the sorting variable (pseudo height of max cooling)
    // for a g point
    std::vector<Real> median_sorting_var;

    for (int jband = 0; jband < nband; ++jband) {
      LOG << "Band " << jband << "\n";
      intVector band_index = find(iband == jband);
      int ibegin = band_index(0);
      int iend   = band_index(end);
      Real min_ref = minval(sorting_variable(range(ibegin,iend)));
      Real max_ref = maxval(sorting_variable(range(ibegin,iend)));
      int ng_band = static_cast<int>((max_ref-min_ref) / max_reflectance_range)+1;
      n_g_points(jband) = ng_band;
      ng += ng_band;

      if (true) {
	// Partition space into equal ranges of reflectance (note that
	// this doesn't account for solar energy in each range)
	Real dref = (max_ref-min_ref) / ng_band + 1.0e-8;
	for (int jg = 0; jg < ng_band; ++jg) {
	  intVector index = find(iband == jband
				 && sorting_variable >= min_ref+jg*dref
				 && sorting_variable <  min_ref+(jg+1)*dref);
	  rank1_per_g_point.push_back(minval(irank(index)));
	  rank2_per_g_point.push_back(maxval(irank(index)));
	  error_per_g_point.push_back(dref);
	  // Make sure that the sorting variables for clouds lie below
	  // those for gases, hence the -2.0
	  median_sorting_var.push_back(-2.0+min_ref+(jg+0.5)*dref);
	  band_num.push_back(jband);
	}
      }
      else {
	// Partition into equal ranges of solar energy
	intVector ireorder(band_index.size());
	ireorder(irank(range(ibegin,iend))-ibegin) = range(ibegin,iend);
	Vector cum_ssi(irank.size());
	cum_ssi = -1.0;
	cum_ssi(ireorder) = cumsum(Vector(ssi(ireorder)));
	Real band_irradiance = sum(ssi(range(ibegin,iend)));
	Real d_irradiance = band_irradiance*(1.0+1.0e-8)/ng_band;
	for (int jg = 0; jg < ng_band; ++jg) {
	  intVector index = find(iband == jband
				 && cum_ssi >= jg*d_irradiance
				 && cum_ssi <  (jg+1)*d_irradiance);
	  rank1_per_g_point.push_back(minval(irank(index)));
	  rank2_per_g_point.push_back(maxval(irank(index)));
	  error_per_g_point.push_back(maxval(sorting_variable(index))
				      -minval(sorting_variable(index)));
	  // Make sure that the sorting variables for clouds lie below
	  // those for gases, hence the -2.0
	  median_sorting_var.push_back(-2.0+mean(sorting_variable(index)));
	  band_num.push_back(jband);
	}
      }
    }

    intVector rank1(&rank1_per_g_point[0], dimensions(ng));
    intVector rank2(&rank2_per_g_point[0], dimensions(ng));
    Vector error(&error_per_g_point[0], dimensions(ng));
    Vector median_sorting_variable(&median_sorting_var[0], dimensions(ng));
    intVector band_number(&band_num[0], dimensions(ng));
    
    SingleGasData cloud_data(cloud_str, n_g_points, band_number,
			     rank1, rank2, error, median_sorting_variable,
			     irank);
    cloud_data.print();
      
    single_gas_data.push_back(cloud_data);
  }


  // Loop over gases
  while (config.read(gas_str, "gases", ngas)) {
    std::string Gas = gas_str;
    std::transform(Gas.begin(), Gas.end(), Gas.begin(), ::toupper);
    LOG << "*** FINDING G POINTS FOR " << Gas << "\n";

    // Used by the shortwave total-transmission averaging method
    Real min_scaling = 1.0, max_scaling = 1.0;
    config.read(min_scaling, gas_str, "min_scaling");
    config.read(max_scaling, gas_str, "max_scaling");
    // We always have to account for some scaling of the optical path
    // due to the variation of solar zenith angle
    min_scaling = std::min(0.5, min_scaling);
    max_scaling = std::max(2.5, max_scaling);

    // READ ORDERING
    std::string reordering_input;
    if (!config.read(reordering_input, gas_str, "reordering_input")) {
      ERROR << "No reordering_input found";
      THROW(PARAMETER_ERROR);
    }
    LOG << "Reading " << reordering_input << "\n";
    DataFile order_file(reordering_input);
    Vector sorting_variable;
    intVector irank, iband;
    order_file.read(irank, "rank");
    order_file.read(iband, "band_number");
    order_file.read(band_bound1, "wavenumber1_band");
    order_file.read(band_bound2, "wavenumber2_band");
    order_file.read(sorting_variable, "sorting_variable");
    nband = band_bound1.size();

    // Read some band-specific configuration for this gas

    // First the number of pieces into which the base g-point in each
    // band should be split after applying the equipartition
    // method. In the near-infrared this is necessary for water
    // vapour.  A value of 1 means do nothing, an integer greater than
    // 1 means split into exactly this number, while a fraction less
    // than 1 means that the number of pieces should be 2 plus the
    // rounded-down result of multiplying this fraction by the
    // existing total number of g-points for this gas and band
    Vector base_split_raw;
    Vector base_split(nband);
    base_split = 1.0;
    if (config.read(base_split_raw, gas_str, "base_split")) {
      int nsize = std::min(nband, base_split_raw.size());
      base_split(range(0,nsize-1)) = base_split_raw(range(0,nsize-1));
      LOG << "Base g-points will be split according to: " << base_split << "\n";
    }

    // The minimum number of g points to use in each band (default 1)
    intVector min_g_points_raw;
    intVector min_g_points(nband);
    min_g_points = 1;
    if (config.read(min_g_points_raw, gas_str, "min_g_points")) {
      int nsize = std::min(nband, min_g_points_raw.size());
      min_g_points(range(0,nsize-1)) = min_g_points_raw(range(0,nsize-1));
    }

    // Set the albedo by band
    band_albedo.resize(nband);
    band_albedo = 0.0;
    intVector iband_no_rayleigh = find(band_bound2 <= max_no_rayleigh_wavenumber);
    band_albedo(iband_no_rayleigh) = reference_albedo;
    max_no_rayleigh_wavenumber = maxval(band_bound2(iband_no_rayleigh));

    // Copy the tolerances to the bandwise vector
    Vector heating_rate_tolerance(nband);
    if (heating_rate_tolerance_in.size() == 1) {
      heating_rate_tolerance = heating_rate_tolerance_in(0);
    }
    else if (heating_rate_tolerance_in.size() == nband) {
      heating_rate_tolerance = heating_rate_tolerance_in;
    }
    else {
      ERROR << "heating_rate_tolerance must have either 1 element, or one per band";
      THROW(PARAMETER_ERROR);
    }

    // "irank" is the rank of each point of the spectrum from the
    // least to the most absorbing.  We want an index that will
    // reorder an array.
    intVector ireorder(irank.size());
    ireorder(irank) = range(0,irank.size()-1);

    sorting_variable = eval(sorting_variable(ireorder));

    // BACKGROUND
  
    Matrix bg_optical_depth;
    // Load background optical depths

    std::string bg_molecules;

    // Wavenumber spacing, cm-1
    Vector d_wavenumber_cm_1;
    // Half-level pressure (Pa) and temperature (K)
    Vector pressure_hl, temperature_hl;
    // Volume mixing ratios of each gas (mol mol-1)
    Matrix vmr_fl;

    bool have_background = false;

    if (config.exist(gas_str + ".background_input")) {
      LOG << "Generating background optical depth\n";

      read_merged_spectrum(config, iprofile, gas_str + ".background_",
			   pressure_hl, temperature_hl,
			   wavenumber_cm_1, d_wavenumber_cm_1,
			   bg_optical_depth, bg_molecules, vmr_fl);
      have_background = true;

      LOG << "  Reordering\n";
      bg_optical_depth = eval(bg_optical_depth(__,ireorder));
    }
    else {
      bg_optical_depth.resize(bg_optical_depth.dimensions());
      bg_optical_depth = 0.0;
    }
  
    // TARGET GAS

    LOG << "Generating target optical depth\n";

    Matrix optical_depth;
    
    read_merged_spectrum(config, iprofile, gas_str + ".",
			 pressure_hl, temperature_hl,
			 wavenumber_cm_1, d_wavenumber_cm_1,
			 optical_depth, molecule, vmr_fl);

    nwav = wavenumber_cm_1.size();

    Vector albedo_orig;
    if (do_sw) {
      albedo_orig.resize(nwav);
      albedo_orig = 0.0;
      albedo_orig(find(wavenumber_cm_1 < max_no_rayleigh_wavenumber)) = reference_albedo;
    }

    LOG << "  Reordering\n";
    optical_depth = eval(optical_depth(__,ireorder));
    wavenumber_cm_1 = eval(wavenumber_cm_1(ireorder));
    d_wavenumber_cm_1 = eval(d_wavenumber_cm_1(ireorder));
 
    int nlay = pressure_hl.size()-1;

    LOG << nlay << " layers\n";
    LOG << nwav << " spectral points\n";

    if (do_sw) {
      albedo = albedo_orig(ireorder);
      // Albedo is constant within a band, so reordering should not do
      // anything - check this!
      if (any(albedo != albedo_orig)) {
	WARNING << "Albedo is not constant within a band";
	ENDWARNING;
      }
      // Check albedo per band
      for (int jband = 0; jband < nband; ++jband) {
	LOG << "  Band " << jband << ": target albedo = " << band_albedo(jband) 
	    << ", mean spectral albedo = "
	    << sum(albedo_orig(find(iband == jband))*ssi(find(iband == jband))) / sum(ssi(find(iband == jband)))
	    << ", mean reordered spectral albedo = "
	    << sum(albedo(find(iband == jband))*ssi(find(iband == jband))) / sum(ssi(find(iband == jband)))
	    << "\n";
      }
    }

    Matrix flux_dn(nlay+1,nwav);
    Matrix flux_up;
    Vector surf_emissivity;

    // For total-transmission method we need to do flux calculations
    // with the target gas optical depth scaled up and down
    Matrix flux_dn_high, flux_dn_low;
    Matrix flux_up_high, flux_up_low;

    if (!do_sw) {

      flux_up.resize(nlay+1,nwav);

      // COMPUTE PLANCK FUNCTION
  
      if (planck_hl.empty()) {

	LOG << "Computing Planck function\n";
  
	planck_hl.resize(nlay+1,nwav);
	planck_function(temperature_hl, wavenumber_cm_1, d_wavenumber_cm_1,
			planck_hl);
	surf_planck.resize(nwav);
	planck_function(temperature_hl(end), wavenumber_cm_1, d_wavenumber_cm_1,
			surf_planck);
      }
      else if (planck_hl.size(0) != nlay+1 || planck_hl.size(1) != nwav) {
	ERROR << "Existing Planck function matrix is wrong size";
	THROW(PARAMETER_ERROR);
      }
      
      // RADIATIVE TRANSFER
    
      surf_emissivity.resize(nwav);
      surf_emissivity = 1.0;
      
      LOG << "Performing longwave radiative transfer\n";
      
      Matrix total_optical_depth = bg_optical_depth + optical_depth;
      radiative_transfer_lw(planck_hl, total_optical_depth, surf_emissivity,
			    surf_planck, flux_dn, flux_up);
      
    }
    else {

      LOG << "Performing shortwave radiative transfer\n";
      // Note that flux_up is not allocated, which ensures that
      // heating rates do not include the upwelling flux contribution

      Matrix total_optical_depth = bg_optical_depth + optical_depth;
      radiative_transfer_direct_sw(cos_sza, ssi,
				   total_optical_depth, flux_dn);

      if (averaging_method == "total-transmission") {
	// We generat two further benchmark spectral radiation fields,
	// one with the optical depths scaled down, the other with
	// them scaled up
	flux_dn_low.resize(nlay+1,nwav);
	flux_dn_high.resize(nlay+1,nwav);
	if (max_no_rayleigh_wavenumber > 0.0) {
	  flux_up_low.resize(nlay+1,nwav);
	  flux_up_high.resize(nlay+1,nwav);
	  total_optical_depth = bg_optical_depth + min_scaling * optical_depth;
	  radiative_transfer_norayleigh_sw(cos_sza, ssi,
					   total_optical_depth, albedo,
					   flux_dn_low, flux_up_low);
	  total_optical_depth = bg_optical_depth + max_scaling * optical_depth;
	  radiative_transfer_norayleigh_sw(cos_sza, ssi,
					   total_optical_depth, albedo,
					   flux_dn_high, flux_up_high);
	}
	else {
	  total_optical_depth = bg_optical_depth + min_scaling * optical_depth;
	  radiative_transfer_direct_sw(cos_sza, ssi,
				       total_optical_depth, flux_dn_low);
	  total_optical_depth = bg_optical_depth + max_scaling * optical_depth;
	  radiative_transfer_direct_sw(cos_sza, ssi,
				       total_optical_depth, flux_dn_high);
	}
      }

    }

    LOG << "Computing heating rate\n";

    Matrix hr(nlay,nwav);
    heating_rate(pressure_hl, flux_dn, flux_up, hr);

    // Save top-of-atmosphere and surface fluxes
    Vector flux_dn_surf, flux_up_toa;
    flux_dn_surf = flux_dn(end,__);

    if (flux_up.empty()) {
      flux_up_toa.resize(flux_dn_surf.size());
      flux_up_toa = 0.0;
    }
    else {
      flux_up_toa  = flux_up(0,__);
    }

    flux_dn.clear();
    flux_up.clear();

    // Additional benchmark values for the shortwave
    // total-transmission averaging method
    Matrix hr_low, hr_high;
    Vector flux_dn_surf_low,  flux_up_toa_low;
    Vector flux_dn_surf_high, flux_up_toa_high;
    if (do_sw && averaging_method == "total-transmission") {
      hr_low.resize(nlay,nwav);
      // We intentionally use the unallocated flux_up array here as
      // the heating rates do not include the upwelling contribution
      heating_rate(pressure_hl, flux_dn_low, flux_up, hr_low);
      flux_dn_surf_low = flux_dn_low(end,__);
      flux_dn_low.clear();
      if (flux_up_low.empty()) {
	flux_up_toa_low.resize(flux_dn_surf_low.size());
	flux_up_toa_low = 0.0;
      }
      else {
	flux_up_toa_low = flux_up_low(0,__);
	flux_up_low.clear();
      }	

      hr_high.resize(nlay,nwav);
      heating_rate(pressure_hl, flux_dn_high, flux_up, hr_high);
      flux_dn_surf_high = flux_dn_high(end,__);
      flux_dn_high.clear();
      if (flux_up_high.empty()) {
	flux_up_toa_high.resize(flux_dn_surf_high.size());
	flux_up_toa_high = 0.0;
      }
      else {
	flux_up_toa_high = flux_up_high(0,__);
	flux_up_high.clear();
      }	
    }

    Vector layer_weight = sqrt(pressure_hl(range(1,end)))-sqrt(pressure_hl(range(0,end-1)));
    Vector pressure_fl = 0.5*(pressure_hl(range(1,end))+pressure_hl(range(0,end-1)));
    
    Real min_pressure = 0.0;
    config.read(min_pressure, "min_pressure");
    layer_weight(find(pressure_fl < min_pressure)) = 0.0;
    layer_weight /= sum(layer_weight);

    // Find "importance" of a particular wavenumber as the RMS heating rate 
    Vector rms_hr = sqrt(sum(spread<1>(layer_weight,nwav)*hr*hr,0));

    LOG << "Finding g points: ";

    std::vector<int> n_g_points_per_band;
    // Lower and upper g
    std::vector<int> rank1_per_g_point, rank2_per_g_point, band_num;
    // Heating rate error, K d-1
    std::vector<Real> error_K_d_1;
    // Median of the sorting variable (pseudo height of max cooling)
    // for a g point
    std::vector<Real> median_sorting_var;

    // G point to which each wavenumber is assigned (unused)
    //intVector g_point(nwav);
    //g_point = -1;

    Matrix metric;
    if (averaging_method == "linear"
	|| averaging_method == "logarithmic"
	|| averaging_method == "total-transmission") {
      metric >>= optical_depth;
      if (averaging_method == "linear") {
	LOG << "linear averaging of optical depth\n";
      }
      else if (averaging_method == "logarithmic") {
	LOG << "logarithmic averaging of optical depth\n";
      }
      else {
	LOG << "Total-transmission averaging of optical depth, scaling "
	    << min_scaling << "-" << max_scaling << "\n";
      }
    }
    else if (averaging_method == "transmission") {
      LOG << "transmission averaging of optical depth\n";
      metric = 1.0-exp(-optical_depth*LW_DIFFUSIVITY);
    }
    else if (averaging_method == "transmission-2") {
      LOG << "2x transmission averaging of optical depth\n";
      metric = 1.0-exp(-optical_depth*LW_DIFFUSIVITY*2.0);
    }
    else if (averaging_method == "square-root") {
      LOG << "square-root averaging of optical depth\n";
      metric = sqrt(optical_depth);
    }
    else {
      ERROR << "Averaging method \"" << averaging_method << "\" not understood";
      THROW(PARAMETER_ERROR);
    }
    
    for (int jband = 0; jband < nband; ++jband) {
      LOG << "Band " << jband << "\n";
      intVector band_index = find(iband == jband);
      int ibegin = band_index(0);
      int iend   = band_index(end);
      //      intVector rank_band = irank(band_index);
      //std::cout << "ibegin = " << ibegin << "  iend = " << iend << std::endl;
      CkdEquipartition Eq;

      if (!do_sw) {
	Eq.init_lw(averaging_method, flux_weight, layer_weight,
		   pressure_hl, surf_emissivity,
		   surf_planck, flux_dn_surf, flux_up_toa, planck_hl,
		   bg_optical_depth, metric, hr, ibegin, iend);
      }
      else {
	Eq.init_sw(averaging_method, flux_weight, layer_weight,
		   cos_sza, pressure_hl, ssi, band_albedo(jband),
		   flux_dn_surf, flux_up_toa,
		   bg_optical_depth, metric, hr, ibegin, iend);
	if (averaging_method == "total-transmission") {
	  Eq.init_sw_extras(flux_dn_surf_low, flux_up_toa_low,
			    flux_dn_surf_high, flux_up_toa_high,
			    min_scaling, max_scaling, 
			    hr_low, hr_high, ibegin, iend);
	}
      }
      Eq.set_partition_max_iterations(max_iterations);
      Eq.set_partition_tolerance(tolerance_tolerance);
      int ng = 10;
#define PARTITION_BY_ERROR 1
#ifdef PARTITION_BY_ERROR
      std::vector<ep_real> bounds, error;
      EpStatus istatus = Eq.equipartition_e(heating_rate_tolerance(jband), 
					    0.0, 1.0, ng, bounds, error);
      if (ng < min_g_points(jband)) {
	ep_print_result(istatus, 1, ng, &bounds[0], &error[0]);
	std::cout << "      computational cost = " << Eq.total_comp_cost << "\n";
	LOG << "  " << ng << " intervals is fewer than minimum of " << min_g_points(jband) << "\n";
	ng = min_g_points(jband);
	bounds.resize(ng+1);
	error.resize(ng);
	// Set initial bounds
	for (int ibound = 0; ibound < ng+1; ++ibound) {
	  bounds[ibound] = sqrt(static_cast<Real>(ibound) / static_cast<Real>(ng));
	}
	istatus = Eq.equipartition_n(ng, &bounds[0], &error[0]);	
      }
#else
      Vector bounds(ng+1), error(ng);
      bounds = sqrt(linspace(0.0, 1.0, ng+1));
      EpStatus istatus = Eq.equipartition_n(ng, &bounds[0], &error[0]);
#endif

      ep_print_result(istatus, 1, ng, &bounds[0], &error[0]);
      std::cout << "      computational cost = " << Eq.total_comp_cost << "\n";

      if (base_split(jband) != 1.0) {
	// Work out how many pieces to split the base interval into
	int nsplit = 1;
	if (base_split(jband) > 1.0) {
	  nsplit = base_split(jband);
	  if (nsplit == 1) {
	    ERROR << "Positive values of base_split must be at least 2";
	    THROW(PARAMETER_ERROR);
	  }
	}
	else {
	  // Always split into at least two
	  nsplit = 2+static_cast<int>(base_split(jband)*ng);
	}

	LOG << "  Splitting base interval into " << nsplit << " pieces\n";
	Real first_bound = bounds[1];
	// First error is now incorrect
	error[0] = -1.0;
	for (int ibnd = 0; ibnd < nsplit-1; ++ibnd) {
	  bounds.insert(bounds.begin()+ibnd+1, first_bound*(ibnd+1)/static_cast<Real>(nsplit));
	  error.insert(error.begin()+ibnd, -1.0);
	}
	ng += nsplit-1;
	ep_print_result(istatus, 1, ng, &bounds[0], &error[0]);
      }

      Vector bnds(&bounds[0], dimensions(bounds.size()));
      Vector err(&error[0], dimensions(error.size()));
      std::cout << "      bounds = " << bnds << "\n";
      std::cout << "      error  = " << err << "\n";
      n_g_points_per_band.push_back(ng);
      for (int ig = 0; ig < ng; ++ig) {
	int ind1 = Eq.lower_index(bounds[ig])+ibegin;
	int ind2 = Eq.upper_index(bounds[ig+1])+ibegin;

	rank1_per_g_point.push_back(ind1);
	rank2_per_g_point.push_back(ind2);
	error_K_d_1.push_back(error[ig]);
	band_num.push_back(jband);
	if (!do_sw) {
	  median_sorting_var.push_back(calc_median_sorting_variable(sorting_variable, surf_planck, ind1, ind2));
	}
	else {
	  median_sorting_var.push_back(calc_median_sorting_variable(sorting_variable, ssi, ind1, ind2));
	}

	//g_point(rank_band(range(ind1,ind2))) = ig;
	// This should probably be ireorder rather than irank, but it's not used anyway
	//g_point(irank(range(ind1,ind2))) = ig;
      }

      if (debug_partition) {
	std::cerr << "debug_partition_" << Gas << "_band = " << jband << "\n";
	Eq.debug_partition = true;
	Eq.set_parallel(false);
	Eq.calc_error_all(ng, &bounds[0], &error[0]);
	Eq.set_parallel(true);
	Eq.debug_partition = false;
      }
    
    } // Loop over bands
    metric.clear();

    // Create intVector pointing to std::vector
    intVector n_g_points(&n_g_points_per_band[0], dimensions(nband));
    int ng = rank1_per_g_point.size();

    intVector rank1(&rank1_per_g_point[0], dimensions(ng));
    intVector rank2(&rank2_per_g_point[0], dimensions(ng));
    Vector error(&error_K_d_1[0], dimensions(ng));
    Vector median_sorting_variable(&median_sorting_var[0], dimensions(ng));
    intVector band_number(&band_num[0], dimensions(ng));

    SingleGasData raw_single_gas_data(gas_str, n_g_points, band_number,
				      rank1, rank2, error, median_sorting_variable,
				      irank);
    raw_single_gas_data.print();

    single_gas_data.push_back(raw_single_gas_data);
    ++ngas;
    LOG << "\n";
  } // Loop over gases

  // If cloud is present then it behaves as a pseudo-gas
  if (config.read(cloud_str, "cloud")) {
    ngas++;
  }

  LOG << "*** COMPUTING SPECTRAL OVERLAP OF GASES\n";
  // "ng" is the total number of g points for all gases
  intVector band_number;
  int ng = overlap_g_points(single_gas_data, band_number);

  // Work out to which multi-gas g point each wavenumber is assigned
  intVector g_point(nwav);
  
  g_point = -1;
  boolVector is_found(nwav);
  for (int ig = 0; ig < ng; ++ig) {
    is_found = true; // Start with each wavenumber being "in" the g
		     // point then progressively eliminate wavenumbers
    for (int igas = 0; igas < ngas; ++igas) {
      is_found.where(  single_gas_data[igas].g_point < single_gas_data[igas].g_min(ig)
		     ||single_gas_data[igas].g_point > single_gas_data[igas].g_max(ig)) = false;
    }
    g_point.where(is_found) = ig;
    if (!any(is_found)) {
      WARNING << "g point " << ig << " occupies none of the spectrum";
      ENDWARNING;
    }
  }
  is_found.clear();

  if (any(g_point == -1)) {
    WARNING << count(g_point == -1) << " wavenumbers are not assigned to a g point";
    ENDWARNING;
  }

  LOG << "\n";

  LOG << "Writing " << output << "\n";

  OutputDataFile file(output);

  // Define dimensions

  file.define_dimension("band", nband);
  if (ng > 0) {
    file.define_dimension("g_point", ng);
  }

  std::string molecule_list;
  for (int igas = 0; igas < ngas; ++igas) {
    file.define_dimension(single_gas_data[igas].molecule + "_g_point", 
			  single_gas_data[igas].rank1.size());
    if (igas == 0) {
      molecule_list = single_gas_data[igas].molecule;
    }
    else {
      molecule_list += " " + single_gas_data[igas].molecule;
    }
  }
  
  if (nwav > 0) {
    file.define_dimension("wavenumber", nwav);
  }

  // Define variables

  file.define_variable("n_gases", INT);
  file.write_long_name("Number of gases treated", "n_gases");
  file.write_comment("The gases are listed in the global attribute \"constituent_id\".", "n_gases");
  
  file.define_variable("wavenumber1_band", FLOAT, "band");
  file.write_long_name("Lower wavenumber bound of band", "wavenumber1_band");
  file.write_units("cm-1", "wavenumber1_band");

  file.define_variable("wavenumber2_band", FLOAT, "band");
  file.write_long_name("Upper wavenumber bound of band", "wavenumber2_band");
  file.write_units("cm-1", "wavenumber2_band");

  file.define_variable("band_number", SHORT, "g_point");
  file.write_long_name("Band number of each g point", "band_number");
  
  if (!ssi.empty()) {
    // This is a shortwave file
    file.define_variable("solar_irradiance", FLOAT, "g_point");
    file.write_long_name("Solar irradiance across each g point", "solar_irradiance");
    file.write_units("W m-2", "solar_irradiance");
  }

  for (int igas = 0; igas < ngas; ++igas) {
    const SingleGasData& this_gas = single_gas_data[igas];
    const std::string& Molecule = this_gas.Molecule;
    const std::string& molecule = this_gas.molecule;
    std::string dim_str = molecule + "_g_point";
    file.define_variable(molecule + "_n_g_points", INT, "band");
    file.write_long_name(std::string("Number of g points for ")
			 + Molecule, molecule + "_n_g_points");

    file.define_variable(molecule + "_band_number", SHORT, dim_str);
    file.write_long_name("Band number of each " + Molecule + " g point",  molecule + "_band_number");
    file.write_comment("This variable indicates the number of the band (0 based) that each g point is in.", molecule + "_band_number");

    file.define_variable(molecule + "_rank1", INT, dim_str);
    file.write_long_name("Rank of first wavenumber for " + Molecule, molecule + "_rank1");

    file.define_variable(molecule + "_rank2", INT, dim_str);
    file.write_long_name("Rank of last wavenumber for " + Molecule, molecule + "_rank2");

    file.define_variable(molecule + "_error", FLOAT, dim_str);
    file.write_long_name("Root-mean-square heating-rate error for " + Molecule, molecule + "_error");
    file.write_units("K d-1", molecule + "_error");

    file.define_variable(molecule + "_sorting_variable", FLOAT, dim_str);
    file.write_long_name(std::string("Median in g-point of variable used to sort ")
			 + Molecule + " spectrum", molecule + "_sorting_variable");

    if (ng > 0) {
      file.define_variable(molecule + "_g_min", INT, "g_point");
      file.write_long_name(std::string("Minimum ") + Molecule 
			   + " g point contributing to merged g points", molecule + "_g_min");
      
      file.define_variable(molecule + "_g_max", INT, "g_point");
      file.write_long_name(std::string("Maximum ") + Molecule 
			   + " g point contributing to merged g points", molecule + "_g_max");
    }
  }

  if (nwav > 0) {
    file.define_variable("wavenumber", DOUBLE, "wavenumber");
    file.write_long_name("Wavenumber", "wavenumber");
    file.write_units("cm-1", "wavenumber");
    file.define_variable("g_point", SHORT, "wavenumber");
    file.write_long_name("G point", "g_point");
    file.deflate_variable("g_point");
    for (int igas = 0; igas < ngas; ++igas) {
      const SingleGasData& this_gas = single_gas_data[igas];
      const std::string& Molecule = this_gas.Molecule;
      const std::string& molecule = this_gas.molecule;
      file.define_variable(molecule + "_g_point", SHORT, "wavenumber");
      file.write_long_name(Molecule + " g point", molecule + "_g_point");
      file.deflate_variable(molecule + "_g_point");
    }
  }

  // Define global variables
  
  std::string title;
  if (ssi.empty()) {
    title = "Definition of the spectral intervals of a longwave CKD model";
  }
  else {
    title = "Definition of the spectral intervals of a shortwave CKD model";
  }
 
  write_standard_attributes(file, title);

  file.write(molecule_list, "constituent_id");

  file.append_history(argc, argv);

  std::string config_str;
  config.read(config_str);  
  file.write(config_str, "config");

  // Write data

  file.end_define_mode();

  file.write(ngas, "n_gases");

  file.write(band_bound1, "wavenumber1_band");
  file.write(band_bound2, "wavenumber2_band");

  file.write(band_number, "band_number");
  if (!ssi.empty()) {
    Vector solar_irradiance(ng);
    solar_irradiance = 0.0;
    for (int ig = 0; ig < ng; ++ig) {
      solar_irradiance(ig) = sum(ssi(find(g_point == ig)));
    }
    int nbad = count(solar_irradiance <= 0.0);
    if (nbad > 0) {
      WARNING << nbad << " shortwave g points have zero solar irradiance";
      ENDWARNING;
    }
    file.write(solar_irradiance, "solar_irradiance");
  }

  for (int igas = 0; igas < ngas; ++igas) {
    const SingleGasData& this_gas = single_gas_data[igas];
    const std::string& molecule = this_gas.molecule;
    file.write(this_gas.n_g_points, molecule + "_n_g_points");
    file.write(this_gas.band_number, molecule + "_band_number");
    file.write(this_gas.rank1, molecule + "_rank1");
    file.write(this_gas.rank2, molecule + "_rank2");
    file.write(this_gas.error, molecule + "_error");
    file.write(this_gas.sorting_variable, molecule + "_sorting_variable");
    if (ng > 0) {
      file.write(this_gas.g_min, molecule + "_g_min");
      file.write(this_gas.g_max, molecule + "_g_max");
    }
  }

  if (nwav > 0) {
    file.write(wavenumber_cm_1, "wavenumber");
    file.write(g_point, "g_point");
    for (int igas = 0; igas < ngas; ++igas) {
      const SingleGasData& this_gas = single_gas_data[igas];
      const std::string& molecule = this_gas.molecule;
      file.write(this_gas.g_point, molecule + "_g_point");
    }
  }

  file.close();

}
