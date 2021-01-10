#include <string>
#include <algorithm>

#include "OutputDataFile.h"
#include "write_order.h"

/// Write a NetCDF file containing the ordering of spectral intervals
/// for a single gas for use in a correlated k-distribution scheme
void
write_order(std::string& file_name,                    ///< Name of NetCDF file to write
	    int argc,                                  ///< Number of command-line args
	    const char** argv,                         ///< Command-line arguments
	    std::string& molecule,                     ///< Formula for gas in lower case
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
	    ) {
  int nband= band_bound1.size();
  int nwav = wavenumber_cm_1.size();
  OutputDataFile file(file_name);

  // Define dimensions

  file.define_dimension("band", nband);
  file.define_dimension("wavenumber", nwav);

  // Define variables

  file.define_variable("wavenumber1_band", FLOAT, "band");
  file.write_long_name("Lower wavenumber bound of band", "wavenumber1_band");
  file.write_units("cm-1", "wavenumber1_band");

  file.define_variable("wavenumber2_band", FLOAT, "band");
  file.write_long_name("Upper wavenumber bound of band", "wavenumber2_band");
  file.write_units("cm-1", "wavenumber2_band");

  file.define_variable("wavenumber", DOUBLE, "wavenumber");
  file.deflate_variable("wavenumber");
  file.write_long_name("Wavenumber", "wavenumber");
  file.write_units("cm-1", "wavenumber");

  file.define_variable("d_wavenumber", FLOAT, "wavenumber");
  file.deflate_variable("d_wavenumber");
  file.write_long_name("Wavenumber interval", "d_wavenumber");
  file.write_units("cm-1", "d_wavenumber");

  file.define_variable("band_number", SHORT, "wavenumber");
  file.deflate_variable("band_number");
  file.write_long_name("Band number", "band_number");
  file.write_comment("This variable indicates the number of the band (0 based) that each wavenumber is in, with -1 indicating a wavenumber not considered.", "band_number");

  file.define_variable("rank", INT, "wavenumber");
  file.deflate_variable("rank");
  file.write_long_name("Rank when reordered", "rank");
  file.write_comment("This variable indicates the place of each wavenumber after reordering, with 0 indicating the least optically thick.\n"
		     "rank(i) provides the rank of wavenumber i.", "rank");

  /*
  file.define_variable("ordered_index", INT, "wavenumber");
  file.deflate_variable("ordered_index");
  file.write_long_name("Index to ordered points", "ordered_index");
  file.write_comment("ordered_index(i) provides the index to the ith least optically thick wavenumber.\n"
		     "It is related to the rank variable via rank(ordered_index(i))=i.", "ordered_index");
  */

  file.define_variable("column_optical_depth", FLOAT, "wavenumber");
  file.deflate_variable("column_optical_depth");
  file.write_long_name("Column optical depth", "column_optical_depth");

  file.define_variable("sorting_variable", FLOAT, "wavenumber");
  file.deflate_variable("sorting_variable");
  file.write_long_name("Variable used to sort spectrum", "sorting_variable");
  file.write_comment("This variable is equal to log(surface pressure) minus log(pressure of peak heating/cooling),\n"
		     "but for column optical depths less than a threshold, set to column optical depth minus the threshold.",
		     "sorting_variable");

  // Define global variables

  if (!molecule.empty()) {
    std::string Molecule = molecule;
    std::transform(Molecule.begin(), Molecule.end(), Molecule.begin(), ::toupper);
    std::string title = "Optimal reordering of the absorption spectrum of ";
    title += Molecule;
    file.write(title, "title");

    file.write(molecule, "molecule");
  }
  else {
    std::string title = "Optimal reordering of the absorption spectrum of a gas";
    file.write(title, "title");
  }
  file.append_history(argc, argv);
  file.write(config_str, "config");


  // Write data

  file.end_define_mode();

  file.write(band_bound1, "wavenumber1_band");
  file.write(band_bound2, "wavenumber2_band");
  file.write(wavenumber_cm_1, "wavenumber");
  file.write(d_wavenumber_cm_1, "d_wavenumber");
  file.write(iband, "band_number");
  file.write(rank, "rank");
  //  file.write(ordered_index, "ordered_index");
  file.write(column_optical_depth, "column_optical_depth");
  file.write(peak_cooling_height, "sorting_variable");

  file.close();

}
