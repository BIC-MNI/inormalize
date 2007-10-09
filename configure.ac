# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([inormalize],[1.0.2], [Alex Zidjenbos <alex@bic.mni.mcgill.ca>])
AC_CONFIG_SRCDIR([inormalize.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
# AC_DISABLE_SHARED

# figure out where perl is
AC_PATH_PROGS(PERL, [perl], no)

# Checks for libraries.  See m4/README.
mni_REQUIRE_VOLUMEIO
mni_REQUIRE_EBTKS

# AC_PROG_LIBTOOL

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT( 
epm-header
nu_correct_norm
normalize_mri
lgmask
headmask
)
