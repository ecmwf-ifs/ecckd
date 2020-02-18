# ---------------------------------------------------------------------------
# COPYRIGHT    : Alessio Bozzo - ECMWF
# ----------------------------------------------------------------------------
# PROJECT      : DORSY
# FILE         : netcdf.m4
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

# AX_CHECK_NETCDF([ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Upon success, sets the variable NETCDF_LDFLAGS.

dnl defines a custom macro
AC_DEFUN([AX_CHECK_NETCDF], [

      dnl provides a framework to handle the --with-{arg} values passed to configure on the command line      
      AC_ARG_WITH([netcdf],
            [AS_HELP_STRING([--with-netcdf=DIR], [use NetCDF Library from directory DIR])],
            netcdf_prefix="$with_netcdf"
            []
            )
      
      AS_IF([test x$netcdf_prefix != x],
            [AS_IF([test -d "$netcdf_prefix/lib"],
                  [NETCDF_LDFLAGS="-L$netcdf_prefix/lib -Wl,-rpath,$netcdf_prefix/lib"
                  NETCDF_CPPFLAGS="-I$netcdf_prefix/include"],
		  [test -d "$netcdf_prefix/lib64"],
                  [NETCDF_LDFLAGS="-L$netcdf_prefix/lib64 -Wl,-rpath,$netcdf_prefix/lib64"
                  NETCDF_CPPFLAGS="-I$netcdf_prefix/include"],
                  [AC_MSG_ERROR([
  -----------------------------------------------------------------------------
     --with-netcdf=$netcdf_prefix is not a valid directory
  -----------------------------------------------------------------------------])])],
      [AC_MSG_WARN([
  -----------------------------------------------------------------------------
   Missing option `--with-netcdf=DIR`. Looking for NetCDF Library
   into Linux default library search paths
  -----------------------------------------------------------------------------])]
           )
     
      LDFLAGS="$NETCDF_LDFLAGS $LDFLAGS"
      CPPFLAGS="$NETCDF_CPPFLAGS $CPPFLAGS"
      ax_have_netcdf=yes
      dnl checks for NETCDF
      AC_SEARCH_LIBS([nc_create], [netcdf], [], AC_MSG_ERROR([Unable to find NetCDF library]))
      

      AS_IF([test "x$ax_have_netcdf" = xyes],
            dnl outputing NetCDF Library
            [AC_SUBST([NETCDF_LDFLAGS])
            AC_SUBST([NETCDF_CPPFLAGS])
            $1],
            [$2])
      ]
)
dnl vim:set softtabstop=4 shiftwidth=4 expandtab:
