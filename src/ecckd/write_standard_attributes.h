// write_standard_attributes.h - Write standard global attributes to a NetCDF file
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

#ifndef WRITE_STANDARD_ATTRIBUTES_H
#define WRITE_STANDARD_ATTRIBUTES_H 1

#include "config.h"
#include "OutputDataFile.h"

inline 
void write_standard_attributes(OutputDataFile& file, const std::string& title) {
  file.write(title, "title");
  file.write(PACKAGE_NAME " gas optics tool", "source");
  file.write(PACKAGE "-" PACKAGE_VERSION, "source_id");
  file.write(PACKAGE_VERSION, "software_version");
}

#endif
