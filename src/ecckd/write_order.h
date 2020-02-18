#include <string>
#include <adept_arrays.h>

/// Write a NetCDF file containing the ordering of spectral intervals
/// for a single gas for use in a correlated k-distribution scheme
void
write_order(std::string& file_name,                    ///< Name of NetCDF file to write
	    int argc,                                  ///< Number of command-line args
	    const char** argv,                         ///< Command-line arguments
	    std::string& gas_name,                     ///< Formula for gas in lower case
	    std::string& config_str,                   ///< Configuration as a single string
	    const adept::Vector& band_bound1,          ///< Lower wavenumber of bands (cm-1)
	    const adept::Vector& band_bound2,          ///< Upper wavenumber of bands (cm-1)
	    const adept::Vector& wavenumber_cm_1,      ///< Wavenumber (cm-1)
	    const adept::Vector& d_wavenumber_cm_1,    ///< Wavenumber interval (cm-1)
	    const adept::intVector& iband,             ///< Band number (0 based)
	    const adept::intVector& rank,              ///< Rank of point (0 based)
	    const adept::intVector& ordered_index,     ///< Index to points in order (0 based)
	    const adept::Vector& column_optical_depth, ///< Column optical depth
	    const adept::Vector& peak_cooling_height   ///< Pseudo height of peak cooling
	    );
