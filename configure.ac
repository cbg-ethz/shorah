dnl   How to cut a new release
dnl   ========================
dnl
dnl   1.) Bootstrap the autotools by running:
dnl
dnl           autoreconf -vif
dnl
dnl   2.) Run the configure script:
dnl
dnl           ./configure
dnl
dnl   3.) Check that the build system and out-of-source builds work fine:
dnl
dnl           make distcheck
dnl
dnl   4.) If the last command succeeded successfully, upload the produced bz2 tarball to the Github releases page:
dnl
dnl           https://github.com/cbg-ethz/shorah/releases


dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[shorah],
	[1.9.99],
	[osvaldo.zagordi@gmail.com])
AC_CONFIG_SRCDIR([src/cpp/dpm_sampler.cpp])
AC_CONFIG_MACRO_DIR([m4])



dnl ================================================================================
dnl Check whether we want to set defaults for CFLAGS, CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ================================================================================
AC_MSG_CHECKING([whether configure should try to set CXXFLAGS/CPPFLAGS/LDFLAGS])
AS_IF([test "x${CFLAGS+set}" = "xset" || test "x${CXXFLAGS+set}" = "xset" || test "x${CPPFLAGS+set}" = "xset" || test "x${LDFLAGS+set}" = "xset"],
	[enable_flags_setting=no],
	[enable_flags_setting=yes]
)
AC_MSG_RESULT([${enable_flags_setting}])

AX_CHECK_ENABLE_DEBUG
AS_IF([test "x${enable_debug}" = "xno"],
	[
		cflags_test="-ffast-math -O3"
		cppflags_test=""

		AC_CANONICAL_HOST
		AS_CASE([${host_os}],
			[darwin*],
			[ldflags_test="-Wl,-dead_strip_dylibs"],
			[linux*],
			[ldflags_test="-Wl,-O1 -Wl,--as-needed"]
		)
	]
)



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AC_PROG_SED
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR

AC_LANG([C])
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

AC_LANG([C++])
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

dnl required by samtools on linux
AX_PTHREAD



dnl ==================================
dnl Set CXXFLAGS, CPPFLAGS and LDFLAGS
dnl ==================================
AS_IF([test "x${enable_flags_setting}" = "xyes" && test "x${enable_debug}" = "xno"], [
	AC_LANG([C])
	AX_APPEND_COMPILE_FLAGS([${cflags_test}], [CFLAGS])
	CFLAGS=$( echo ${CFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

	AC_LANG([C++])
	AX_APPEND_COMPILE_FLAGS([${cflags_test}], [CXXFLAGS])
	CXXFLAGS=$( echo ${CXXFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

	AX_APPEND_COMPILE_FLAGS([${cppflags_test}], [CPPFLAGS])
	CPPFLAGS=$( echo ${CPPFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )

	AX_APPEND_LINK_FLAGS([${ldflags_test}], [LDFLAGS])
	LDFLAGS=$( echo ${LDFLAGS} | $SED -e 's/^ *//' -e 's/ *$//' )
])

AC_ARG_ENABLE([popcnt],
	AS_HELP_STRING([--disable-popcnt], [Disable using POPCNT for Hamming distances]))

dnl Check for popcount intrinsics
AS_IF([test "x$enable_popcnt" != "xno"], [
	AC_LANG([C++])
	AX_APPEND_COMPILE_FLAGS([-mpopcnt], [CXXFLAGS])

	AC_MSG_CHECKING([whether configure can enable POPCNT for Hamming distances])
	AC_LINK_IFELSE([
		AC_LANG_SOURCE([[
			#include <nmmintrin.h>
			int main() { return _mm_popcnt_u64(0); }
		]])
	], [
		dnl POPCNT works
		AC_DEFINE([HAVE_POPCNT], [1], [Define if popcnt intrinsics are available])
		AC_MSG_RESULT([yes])
	], [
		dnl POPCNT is broken
		AC_MSG_RESULT([no])
	])
])



dnl =================
dnl Find zlib and gsl
dnl =================
PKG_CHECK_MODULES([ZLIB], [zlib])
PKG_CHECK_MODULES([GSL], [gsl])



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.15 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])
AM_SILENT_RULES([yes])



dnl ======================================
dnl Initialise Perl and Python directories
dnl ======================================
AM_PATH_PYTHON([3.5])



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	${PACKAGE_NAME} ${VERSION}

	CC : ...................... ${CC} (${ax_cv_c_compiler_vendor}, ${ax_cv_c_compiler_version})
	CFLAGS : .................. ${CFLAGS}

	CXX : ..................... ${CXX} (${ax_cv_cxx_compiler_vendor}, ${ax_cv_cxx_compiler_version})
	CXXFLAGS : ................ ${CXXFLAGS}

	CPPFLAGS : ................ ${CPPFLAGS}
	LDFLAGS : ................. ${LDFLAGS}

	PYTHON : .................. ${PYTHON}
	Python version : .......... ${PYTHON_VERSION}
])
