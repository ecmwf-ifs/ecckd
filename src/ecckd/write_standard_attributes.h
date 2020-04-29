#ifndef WRITE_STANDARD_ATTRIBUTES_H
#define WRITE_STANDARD_ATTRIBUTES_H 1

#include "config.h"
#include "OutputDataFile.h"

void write_standard_attributes(OutputDataFile& file, const std::string& title) {
  file.write(title, "title");
  file.write(PACKAGE_NAME " gas optics tool", "source");
  file.write(PACKAGE "-" PACKAGE_VERSION, "source_id");
  file.write(PACKAGE_VERSION, "software_version");
}

#endif
