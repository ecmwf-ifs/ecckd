// merge_spectra.cpp - Merge spectra from multiple gases
//
// Copyright (C) 2019- ECMWF.
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

#include "read_merged_spectrum.h"
#include "DataFile.h"
#include "OutputDataFile.h"

int
main(int argc, const char* argv[])
{
  using namespace adept;

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

  // Names of input and output files
  std::string input, output;

  // File containing the optical depth of the target gas
  if (!config.read(input, "input")) {
    ERROR << "\"input\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  if (!config.read(output, "output")) {
    ERROR << "\"output\" file not specified";
    THROW(PARAMETER_ERROR);
  }

  int ngas, ncol;
  std::string molecules;

  Matrix optical_depth;

  // Wavenumber and wavenumber spacing, cm-1
  Vector wavenumber_cm_1, d_wavenumber_cm_1;
  // Half-level pressure (Pa) and temperature (K)
  Vector pressure_hl, temperature_hl;
  // Volume mixing ratio of each gas (mol mol-1)
  Matrix vmr_fl;

  int icol = 0;

  LOG << "Merging profile " << icol << "\n";

  read_merged_spectrum(config, icol, "",
		       pressure_hl, temperature_hl,
		       wavenumber_cm_1, d_wavenumber_cm_1,
		       optical_depth, molecules, vmr_fl, &ngas, &ncol);

  int nlay = optical_depth.size(0);
  int nwav = optical_depth.size(1);

  LOG << "Writing " << output << "\n";

  OutputDataFile file(output);

  // Define dimensions

  file.define_dimension("column", ncol);
  file.define_dimension("level", nlay);
  file.define_dimension("half_level", nlay+1);
  file.define_dimension("wavenumber", nwav);

  // Define variables

  file.define_variable("pressure_hl", FLOAT, "column", "half_level");
  file.write_long_name("Pressure at half levels", "pressure_hl");
  file.write_units("Pa", "pressure_hl");

  file.define_variable("temperature_hl", FLOAT, "column", "half_level");
  file.write_long_name("Temperature at half levels", "temperature_hl");
  file.write_units("K", "temperature_hl");

  file.define_variable("wavenumber", DOUBLE, "wavenumber");
  file.deflate_variable("wavenumber");
  file.write_long_name("Wavenumber", "wavenumber");
  file.write_units("cm-1", "wavenumber");

  file.define_variable("optical_depth", FLOAT, "column", "level", "wavenumber");
  file.deflate_variable("optical_depth");
  file.set_chunking({1,nlay,nwav}, "optical_depth");
  file.write_long_name("Layer optical depth", "optical_depth");

  // Define global variables

  std::string Molecules = molecules;
  std::transform(Molecules.begin(), Molecules.end(), Molecules.begin(), ::toupper);
  std::string title = "Merged spectral optical depth profiles of ";
  for (int ic = 0; ic < Molecules.size(); ++ic) {
    if (Molecules[ic] == ',') {
      title += ", ";
    }
    else {
      title += Molecules[ic];
    }
  }

  file.write(title, "title");

  file.write(std::string("hybrid:") + molecules, "molecule");
  file.append_history(argc, argv);

  std::string config_str;
  config.read(config_str);  
  file.write(config_str, "config");

  // Write data for first profile

  file.end_define_mode();

  file.write(pressure_hl, "pressure_hl", icol);
  file.write(temperature_hl, "temperature_hl", icol);
  file.write(wavenumber_cm_1, "wavenumber");
  file.write(optical_depth, "optical_depth", icol);

  for (icol = 1; icol < ncol; ++icol) {
    LOG << "Merging profile " << icol << "\n";
    read_merged_spectrum(config, icol, "",
			 pressure_hl, temperature_hl,
			 wavenumber_cm_1, d_wavenumber_cm_1,
			 optical_depth, molecules, vmr_fl);
    file.write(pressure_hl, "pressure_hl", icol);
    file.write(temperature_hl, "temperature_hl", icol);
    file.write(optical_depth, "optical_depth", icol);
  }
  file.close();
}
