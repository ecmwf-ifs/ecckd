/// @file      Error.h
/// @brief     Handling of fatal exceptions
/// @author    Robin J. Hogan
/// @copyright 2006-2015 University of Reading
/// @copyright (C) Copyright 2015- ECMWF
/// @license   This software is licensed under the terms of the Apache Licence Version 2.0
///            which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
///            In applying this licence, ECMWF does not waive the privileges and immunities
///            granted to it by virtue of its status as an intergovernmental organisation
///            nor does it submit to any jurisdiction.

#ifndef Error_H
#define Error_H 1

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <typeinfo>

#include "EsaExitCodes.h"

// If this is true then THROW will throw an exception that is not
// intended to be caught in order that the stack trace is preserved
extern bool trace_exceptions_;

class Exception {};
//class NonFatalException: public Exception {};
//class ArrayOverflow: public NonFatalException {};
//class ExpressionSizeMismatch: public NonFatalException {};
class FatalException: public Exception {};
//class MemoryAllocationFailure: public FatalException {};
class SegmentationFault: public FatalException {};
class FloatingPointException: public FatalException {};
class Interrupt: public FatalException {};

// This is intended not to be caught: it is thrown when
// trace_exceptions_==true if not caught then the debugger can figure
// out where the error was thrown
class DebugPreserveStackTrace {};
 
inline void set_trace_exceptions(bool b) { trace_exceptions_ = b; }

extern int fatal_exception_status_;

std::string stack_trace();
std::string demangle_cpp_name(const std::string& mangled_name);

//template<class Type> 
//inline
//std::string type_string() { return demangle_cpp_name(typeid(Type).name()); }
#define type_string(Type) demangle_cpp_name(typeid(Type).name())

void install_interrupt_handler();
void install_segmentation_fault_handler();
void install_floating_point_exception_handler();

void handle_interrupt(int code);
void handle_segmentation_fault(int code);
void handle_floating_point_exception(int code);

void check_for_fatal_exception();

void init_error_msg_map();
std::string error_code_msg(int code);


// Formerly the logging stuff was in this file
#include "Logging.h"


#endif
