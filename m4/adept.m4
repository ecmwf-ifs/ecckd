# ---------------------------------------------------------------------------
# COPYRIGHT    : Alessio Bozzo - ECMWF
# ----------------------------------------------------------------------------
# PROJECT      : DORSY
# FILE         : adept.m4
# VERSION      : $Revision: $
# DATE         : $Date: $
# COMPONENT    : ACM-CAP
# TYPE         : m4 file
# IDENTIFIER   : $Id: $
# ----------------------------------------------------------------------------
#
# This file contains a macro processor (m4 file) for configure.ac
#
# ----------------------------------------------------------------------------
#
# LICENSE
# ----------------------------------------------------------------------------
#
# You may obtain a copy of the License distributed with this s/w package, 
# software distributed under the License is distributed on an “AS IS” BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See
# the License for the specific language governing permissions and limitations 
# under the License.
#
# ----------------------------------------------------------------------------

# AX_CHECK_ADEPT([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Upon success, sets the variable ADEPT_LDFLAGS.

dnl defines a custom macro
AC_DEFUN([AX_CHECK_ADEPT], [

      dnl provides a framework to handle the --with-{arg} values passed to configure on the command line      
      AC_ARG_WITH([adept],
            [AS_HELP_STRING([--with-adept=DIR], [use Adept Library from directory DIR])],
            adept_prefix="$with_adept"
            []
            )
      
      AS_IF([test x$adept_prefix != x],
            [AS_IF([test -d "$adept_prefix/lib"],
                  [ADEPT_LDFLAGS="-L$adept_prefix/lib -Wl,-rpath,$adept_prefix/lib"
                  ADEPT_CPPFLAGS="-I$adept_prefix/include"],
		  [test -d "$adept_prefix/lib64"],
                  [ADEPT_LDFLAGS="-L$adept_prefix/lib64 -Wl,-rpath,$adept_prefix/lib64"
                  ADEPT_CPPFLAGS="-I$adept_prefix/include"],
                  [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
     --with-adept=$adept_prefix is not a valid directory
  -----------------------------------------------------------------------------])])],
      [AC_MSG_WARN([
  -----------------------------------------------------------------------------
   Missing option `--with-adept=DIR`. Looking for Adept Library
   into Linux default library search paths
  -----------------------------------------------------------------------------])]
           )
     
      LDFLAGS="$ADEPT_LDFLAGS $LDFLAGS"
      CPPFLAGS="$ADEPT_CPPFLAGS $CPPFLAGS"
      ax_have_adept=yes
      dnl checks for ADEPT
      AC_MSG_CHECKING([for Adept >= 2.0.4: including adept_arrays.h and linking via -ladept])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <adept_arrays.h>
      #include <string>
      #if ADEPT_VERSION < 20004
      #error "Adept version >= 2.0.4 required"
      #endif],[std::string test = adept::compiler_version()])],AC_MSG_RESULT([yes]),AC_MSG_RESULT([no])
      AC_MSG_ERROR([Unable to find Adept library version >= 2.0.4]))

      AS_IF([test "x$ax_have_adept" = xyes],
            dnl outputing Adept Library
            [AC_SUBST([ADEPT_LDFLAGS])
            AC_SUBST([ADEPT_CPPFLAGS])
            $1],
            [$2])
      ]
)
dnl vim:set softtabstop=4 shiftwidth=4 expandtab:
