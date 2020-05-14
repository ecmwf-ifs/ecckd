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
// wavenumbers selected via "index", using the Planck function to
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
// wavenumbers selected via "index", using the spectral solar
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

  }
  void init_sw(std::string am, Real fw, const Vector& lw, Real cs,
	       const Vector& prhl, const Vector& si, 
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
      Vector optical_depth_fit = fit_optical_depth_sw(averaging_method,
						      ibound1, ibound2,
						      ssi.soft_link(),
						      metric.soft_link());
	  
      return calc_cost_function_sw(cos_sza,
				   pressure_hl.soft_link(),
				   ssi.soft_link()(range(ibound1,ibound2)),
				   bg_optical_depth.soft_link()(__,range(ibound1,ibound2)),
				   optical_depth_fit.soft_link(),
				   flux_dn_surf.soft_link()(range(ibound1,ibound2)),
				   flux_up_toa.soft_link()(range(ibound1,ibound2)),
				   hr.soft_link()(__,range(ibound1,ibound2)),
				   flux_weight, layer_weight.soft_link());
    }
  }

  std::string averaging_method;
  Real flux_weight;
  Real cos_sza;
  Vector layer_weight;
  Vector pressure_hl;
  Vector ssi;
  Vector surf_emissivity, surf_planck, flux_dn_surf, flux_up_toa;
  Matrix planck_hl, bg_optical_depth, metric, hr;
  int npoints;
  ep_real total_comp_cost;
  bool do_sw;
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
  Vector ssi;

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

  Real heating_rate_tolerance; // K d-1
  if (!config.read(heating_rate_tolerance, "heating_rate_tolerance")) {
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

  int ngas = 0;
  int nband = 0;
  std::string gas_str;

  Matrix planck_hl;
  Vector surf_planck;
  Vector band_bound1, band_bound2;

  std::vector<SingleGasData> single_gas_data;

  int nwav = 0; // Number of wavenumbers
  // Wavenumber, cm-1
  Vector wavenumber_cm_1;

  // Loop over gases
  while (config.read(gas_str, "gases", ngas)) {
    std::string Gas = gas_str;
    std::transform(Gas.begin(), Gas.end(), Gas.begin(), ::toupper);
    LOG << "*** FINDING G POINTS FOR " << Gas << "\n";

    // READ ORDERING
    std::string reordering_input;
    if (!config.read(reordering_input, gas_str, "reordering_input")) {
      ERROR << "No reordering_input found\n";
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

    LOG << "  Reordering\n";
    optical_depth = eval(optical_depth(__,ireorder));
    wavenumber_cm_1 = eval(wavenumber_cm_1(ireorder));
    d_wavenumber_cm_1 = eval(d_wavenumber_cm_1(ireorder));
 
    int nlay = pressure_hl.size()-1;
    nwav = wavenumber_cm_1.size();

    LOG << nlay << " layers\n";
    LOG << nwav << " spectral points\n";

    Matrix flux_dn(nlay+1,nwav), flux_up(nlay+1,nwav);
    Vector surf_emissivity;

    if (!do_sw) {

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
      
      Matrix total_optical_depth = bg_optical_depth + optical_depth;
      radiative_transfer_direct_sw(cos_sza, ssi,
				   total_optical_depth, flux_dn);
      flux_up = 0.0;

    }

    LOG << "Computing heating rate\n";

    Matrix hr(nlay,nwav);
    heating_rate(pressure_hl, flux_dn, flux_up, hr);

    // Save top-of-atmosphere and surface fluxes
    Vector flux_dn_surf, flux_up_toa;
    flux_dn_surf = flux_dn(end,__);
    flux_up_toa  = flux_up(0,__);

    flux_dn.clear();
    flux_up.clear();

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

    // G point to which each wavenumber is assigned
    intVector g_point(nwav);
    g_point = -1;

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
	LOG << "Total-transmission averaging of optical depth\n";
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
		   cos_sza, pressure_hl, ssi,
		   flux_dn_surf, flux_up_toa,
		   bg_optical_depth, metric, hr, ibegin, iend);
      }
      Eq.set_partition_max_iterations(max_iterations);
      Eq.set_partition_tolerance(tolerance_tolerance);
      int ng = 10;
#define PARTITION_BY_ERROR 1
#ifdef PARTITION_BY_ERROR
      std::vector<ep_real> bounds, error;
      EpStatus istatus = Eq.equipartition_e(heating_rate_tolerance, 
					    0.0, 1.0, ng, bounds, error);
#else
      Vector bounds(ng+1), error(ng);
      bounds = sqrt(linspace(0.0, 1.0, ng+1));
      EpStatus istatus = Eq.equipartition_n(ng, &bounds[0], &error[0]);
#endif

      ep_print_result(istatus, 1, ng, &bounds[0], &error[0]);
      std::cout << "      computational cost = " << Eq.total_comp_cost << "\n";
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
	g_point(irank(range(ind1,ind2))) = ig;
      }
    
    }
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

  file.define_variable("band_number", FLOAT, "g_point");
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

    file.define_variable(molecule + "_band_number", INT, dim_str);
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
    file.write(g_point, "g_point");
    for (int igas = 0; igas < ngas; ++igas) {
      const SingleGasData& this_gas = single_gas_data[igas];
      const std::string& molecule = this_gas.molecule;
      file.write(this_gas.g_point, molecule + "_g_point");
    }
  }

  file.close();

}
