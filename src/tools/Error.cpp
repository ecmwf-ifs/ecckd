/// @file      Error.cpp
/// @brief     Implements handling of fatal exceptions
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#include <iostream>
#include <sstream>
#include <string>
#include <csignal>
#include <cstdlib>
#include <map>

#include <execinfo.h>
#include <cxxabi.h>

#include "Error.h"
#include "EsaExitCodes.h"

// Global variables
int fatal_exception_status_ = 0;
bool trace_exceptions_ = false;
std::map<int,std::string> err_msg_map_;

std::string
demangle_cpp_name(const std::string& mangled_name)
{
  size_t len = 1000;
  char* output = reinterpret_cast<char*>(std::malloc(len*sizeof(char)));
  int status = -999;
  abi::__cxa_demangle(mangled_name.c_str(), output, &len, &status);
  std::string ans;
  if (status == 0) {
    ans = output;
  }
  else if (status == -1) {
    ans = "[Memory allocation error]";
  }
  else if (status == -2) {
    ans = "[Invalid mangled name]";
  }
  else {
    ans = "[An error occurred]";
  }
  if (output) {
    std::free(output);
  }
  return ans;  
}

std::string
stack_trace()
{
  std::stringstream result;
  result << "Stack trace:\n";
  const size_t kMaxDepth = 100;
  void* stackAddrs[kMaxDepth];
  size_t stackDepth;
  char** stackStrings;
  stackDepth = backtrace(stackAddrs, kMaxDepth);
  stackStrings = backtrace_symbols(stackAddrs, stackDepth);
  for( size_t i = 1; i < stackDepth; ++i ) {
    const size_t kMaxNameLen = 4096;
    char function[kMaxNameLen];
    std::string stackString(stackStrings[i]);
    size_t start = stackString.find( '(' );
    size_t end = stackString.find( '+' );
    if (std::string::npos == start || std::string::npos == end) {
      result << stackString;
    }
    int status = 1;
    ++start; // puts us past the '('
    size_t n = end - start;
    size_t len;
    abi::__cxa_demangle(stackString.substr(start, n).c_str(),
			function, &len, &status);
    result << "#" << i << " ";
    if (status == 0) {
      // demanging successful
      result << function << "\n";
    }
    else {
      result << stackString.substr(start, n).c_str() << "()\n";
    }
  }
  result << "(For file names and line numbers, run program in gdb)\n";
  return result.str();
}

void
handle_segmentation_fault(int code)
{
  // Turn off handler to avoid possible infinite recursion
  std::signal(SIGSEGV, SIG_DFL);
  // Report problem
  ERROR << "SEGMENTATION FAULT DETECTED\n";// STAR_LINE;
  ERROR << stack_trace();
  fatal_exception_status_ = SIGSEGV;
}

void
install_segmentation_fault_handler()
{
  std::signal(SIGSEGV, handle_segmentation_fault);
}

void
init_error_msg_map()
{
  err_msg_map_[PREMATURE_TERMINATION]=PREMATURE_TERMINATION_STR;
  err_msg_map_[OUT_OF_MEMORY]=OUT_OF_MEMORY_STR;
  err_msg_map_[UNEXPECTED_EXCEPTION]=UNEXPECTED_EXCEPTION_STR;
  err_msg_map_[DISK_FULL]=DISK_FULL_STR;
  err_msg_map_[XML_ERROR]=XML_ERROR_STR;
  err_msg_map_[XML_WARNING]=XML_WARNING_STR;
  err_msg_map_[MISSING_MANDATORY_FILE]=MISSING_MANDATORY_FILE_STR;
  err_msg_map_[NO_MATCHES_FOR_PATTERN]=NO_MATCHES_FOR_PATTERN_STR;
  err_msg_map_[DUPLICATE_FILE]=DUPLICATE_FILE_STR;
  err_msg_map_[INCONSISTENT_SPACECRAFT_ID]=INCONSISTENT_SPACECRAFT_ID_STR;
  err_msg_map_[UNKNOWN_PROCESSING_TYPE]=UNKNOWN_PROCESSING_TYPE_STR;
  err_msg_map_[NOT_ENOUGH_INPUT_ERROR]=NOT_ENOUGH_INPUT_ERROR_STR;
  err_msg_map_[NOT_ENOUGH_INPUT_WARNING]=NOT_ENOUGH_INPUT_WARNING_STR;
  err_msg_map_[CRITICAL_FRAME_FAILURE]=CRITICAL_FRAME_FAILURE_STR;
  err_msg_map_[RECOVERABLE_FRAME_FAILURE]=RECOVERABLE_FRAME_FAILURE_STR;
  err_msg_map_[PRODUCT_FORMAT_ERROR]=PRODUCT_FORMAT_ERROR_STR;
  err_msg_map_[PRODUCT_FORMAT_WARNING]=PRODUCT_FORMAT_WARNING_STR;
  err_msg_map_[NO_PRODUCT_FOUND_ERROR]=NO_PRODUCT_FOUND_ERROR_STR;
  err_msg_map_[NO_PRODUCT_FOUND_WARNING]=NO_PRODUCT_FOUND_WARNING_STR;
  err_msg_map_[CANNOT_OPEN_MANDATORY_FILE]=CANNOT_OPEN_MANDATORY_FILE_STR;
  err_msg_map_[CANNOT_OPEN_OPTIONAL_FILE]=CANNOT_OPEN_OPTIONAL_FILE_STR;
  err_msg_map_[NO_AUX_FILE_ERROR]=NO_AUX_FILE_ERROR_STR;
  err_msg_map_[NO_AUX_FILE_WARNING]=NO_AUX_FILE_WARNING_STR;
  err_msg_map_[BAD_MANDATORY_AUX_FORMAT]=BAD_MANDATORY_AUX_FORMAT_STR;
  err_msg_map_[BAD_OPTIONAL_AUX_FORMAT]=BAD_OPTIONAL_AUX_FORMAT_STR;
  err_msg_map_[NO_PRODUCT_MODEL]=NO_PRODUCT_MODEL_STR;
  err_msg_map_[READ_ERROR_PRODUCT_MODEL]=READ_ERROR_PRODUCT_MODEL_STR;
  err_msg_map_[FINAL_PRODUCT_CREATION_ERROR]=FINAL_PRODUCT_CREATION_ERROR_STR;
  err_msg_map_[TEMP_FILE_CREATION_ERROR]=TEMP_FILE_CREATION_ERROR_STR;
  err_msg_map_[WRITE_ERROR]=WRITE_ERROR_STR;
  err_msg_map_[PARAMETER_ERROR]=PARAMETER_ERROR_STR;
  err_msg_map_[PARAMETER_WARNING]=PARAMETER_WARNING_STR;
  err_msg_map_[PROCESSING_ERROR]=PROCESSING_ERROR_STR;
  err_msg_map_[PROCESSING_WARNING]=PROCESSING_WARNING_STR;
}

std::string
error_code_msg(int code)
{
  return err_msg_map_[code];
}


void
handle_floating_point_exception(int code)
{
  // Turn off handler to avoid possible infinite recursion
  std::signal(SIGFPE, SIG_DFL);
  // Report problem
  ERROR << "FLOATING POINT EXCEPTION DETECTED\n";// STAR_LINE;
  //  ERROR << stack_trace();
  //  fatal_exception_status_ = SIGFPE;
  THROW(UNEXPECTED_EXCEPTION);
}

void
install_floating_point_exception_handler()
{
  std::signal(SIGFPE, handle_floating_point_exception);
}


void
handle_interrupt(int code)
{
  // Turn off handler to avoid possible infinite recursion
  std::signal(SIGHUP, SIG_DFL);
  std::signal(SIGINT, SIG_DFL);
  std::signal(SIGTERM, SIG_DFL);
  // Report problem
  ERROR << "INTERRUPT DETECTED (CODE " << code << ")\n";// STAR_LINE;
  //  ERROR << stack_trace();

  fatal_exception_status_ = code;
}

void
install_interrupt_handler()
{
  std::signal(SIGHUP, handle_interrupt);
  std::signal(SIGINT, handle_interrupt);
  std::signal(SIGTERM, handle_interrupt);
}

void
check_for_fatal_exception()
{
  if (fatal_exception_status_ == 0) {
    return;
  }
  else if (fatal_exception_status_ == SIGSEGV) {
    throw SegmentationFault();
  }
  else if (fatal_exception_status_ == SIGFPE) {
    throw FloatingPointException();
  }
  else {
    throw Interrupt();
  }
}

