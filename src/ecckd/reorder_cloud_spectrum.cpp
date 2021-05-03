#include <vector>
#include <string>
#include <algorithm>

#include "DataFile.h"
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
  std::string wavenumber_input, input, output;

  // Index to cloud size bin (0-based)
  int isize;

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

  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(isize, "isize")) {
    ERROR << "\"isize\" not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(wavenumber_input, "wavenumber_input")) {
    ERROR << "\"wavenumber_input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  // READ WAVENUMBER GRID
  LOG << "Reading " << wavenumber_input << "\n";

  // Wavenumber and wavenumber spacing, cm-1
  Vector wavenumber_cm_1, d_wavenumber_cm_1;
  {
    DataFile wn_file(wavenumber_input);
    wn_file.read(wavenumber_cm_1, "wavenumber");
  }
  int nwav = wavenumber_cm_1.size();
  d_wavenumber_cm_1.resize(nwav); 
  d_wavenumber_cm_1(range(1,end-1))
    = 0.5 * (wavenumber_cm_1(range(2,end))
	     -wavenumber_cm_1(range(0,end-2)));
  d_wavenumber_cm_1(0) = 0.5*d_wavenumber_cm_1(1);
  d_wavenumber_cm_1(end) = 0.5*d_wavenumber_cm_1(end-1);


  // READ CLOUD PROPERTIES

  LOG << "Reading " << input << "\n";

  Vector ssa, asymmetry, cloud_wavenumber_cm_1;
  {
    DataFile cloud_file(input);
    cloud_file.read(cloud_wavenumber_cm_1, "wavenumber");
    Matrix ssa_mat, asymmetry_mat;
    cloud_file.read(ssa_mat, "single_scattering_albedo");
    cloud_file.read(asymmetry_mat, "asymmetry_factor");
    // Extract the requested size bin
    ssa = ssa_mat[isize];
    asymmetry = asymmetry_mat[isize];
  }

  // Compute absorptance in the large optical depth limit and
  // interpolate on to full wavenumber grid
  Vector abs_inf;
  {
    // First delta-Eddington scaling
    Vector f = asymmetry * asymmetry;
    Vector asymmetry_de = 1.0 / (1.0 + asymmetry);
    Vector ssa_de = ssa * (1.0 - f) / (1.0 - ssa*f);
    Vector abs_inf_c = sqrt((1.0 - ssa_de) / (1.0 - ssa_de*asymmetry_de));
    abs_inf_c = 1.0-(1.0 - abs_inf_c) / (1.0 + abs_inf_c);

    abs_inf = interp(cloud_wavenumber_cm_1, abs_inf_c, wavenumber_cm_1);
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

  std::vector<int> g_index(nwav);
  for (int jwav = 0; jwav < nwav; ++jwav) {
    g_index[jwav] = jwav;
  }
  MyCompare my_compare(abs_inf);

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

  write_order(output, argc, argv, std::string("cloud"), config_str,
	      band_bound_clamp1, band_bound_clamp2,
	      wavenumber_cm_1, d_wavenumber_cm_1,
	      iband, rank, ordered_index, Vector(), abs_inf);

}
