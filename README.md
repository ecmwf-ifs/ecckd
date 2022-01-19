# The ECMWF gas optics tool: ecCKD

This document was last updated 19 January 2022

Author: Robin Hogan <r.j.hogan@ecmwf.int>

For more complete information about compilation and usage of ecCKD,
 please see the documentation on the [ecCKD web
 site](https://confluence.ecmwf.int/x/xwu0dw).

## INTRODUCTION

The ecCKD package can generate fast and accurate gas optics models
suitable for use atmospheric radiation schemes such as
[ecRad](http://confluence.ecmwf.int/display/ECRAD), making use of the
correlated k-distribution (CKD) method, and optionally the
full-spectrum correlated-k (FSCK) method.

## PACKAGE OVERVIEW

The subdirectories are as follows:

- `doc` - Latex documentation

- `src` - C++ source code

- `test` - Scripts for running the various components of ecCKD

- `plot` - Matlab functions for plotting some of the results

- `data` - Cloud scattering file 

- `m4` - Macros supporting the autotools build system

## COMPILING AND RUNNNING ECCKD

For a full description of how to compile and run ecCKD, please consult
the User Guide on the [ecCKD web
site](https://confluence.ecmwf.int/x/xwu0dw).  On a Unix/Linux system
with Latex installed you should be able to recreate the documentation
by typing `make documentation` in the `doc` directory.

## LICENCE

(C) Copyright 2019- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.

Several files in the `src/include` and `src/tools` directories also
contain source files whos copyright is either solely or jointly held
with the University of Reading, but all files are released under the
terms of the same Apache License.

## CONTACT

Please email Robin Hogan <r.j.hogan@ecmwf.int> with any queries or bug
fixes, but note that ECMWF does not commit to providing support to
users of this software.