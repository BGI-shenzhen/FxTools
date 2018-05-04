#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT(FxTools, 0.16, hewm2008@qq.com)
AM_INIT_AUTOMAKE(BamDeal,0.16)
AC_CONFIG_SRCDIR([FxTools.cpp])
AC_CONFIG_HEADER([config.h])
# Checks for programs
AC_PROG_AWK
AC_PROG_CXX([g++])
AC_PROG_CC([gcc])
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main],,AC_MSG_ERROR([You need zlib >= 1.2.3 to build $PACKAGE]))
AC_CHECK_LIB([m], [main],,AC_MSG_ERROR([You need libm to build $PACKAGE]))
AC_CHECK_LIB([ncurses], [main],,AC_MSG_ERROR([You need Pre-install Library ncurses to build $PACKAGE]))
AC_CHECK_LIB([pthread], [main],,AC_MSG_ERROR([You need Pre-install Library pthread to build $PACKAGE]))
AC_CHECK_LIB([boost_system], [main],,AC_MSG_ERROR([You need Pre-install Library boost pthread >= 1.5 to build $PACKAGE]))

AC_CHECK_LIB([bz2],[main])
AC_CHECK_LIB([lzma],[main])
# FIXME: Replace `main' with a function in `-lhts':
AC_CHECK_LIB([hts], [main],,AC_MSG_ERROR([You need Pre-install Library libhts samtools/htslib >= 1.5 to build $PACKAGE]))



# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([zlib.h],[],[AC_MSG_ERROR([You need include zlib.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([ncurses.h],[],[AC_MSG_ERROR([You need include ncurses.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([htslib/sam.h],[],[AC_MSG_ERROR([You need include htslib/sam.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([htslib/kstring.h],[],[AC_MSG_ERROR([You need include htslib/kstring.h to build $PACKAGE])],[])

AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/thread/thread.hpp],[],[AC_MSG_ERROR([You need include Boost to build $PACKAGE])],[])
AC_CHECK_HEADERS([boost/bind.hpp],[],[AC_MSG_ERROR([You need include Boost to build $PACKAGE])],[])
AC_CHECK_HEADERS([boost/thread/mutex.hpp],[],[AC_MSG_ERROR([You need include Boost to build $PACKAGE])],[])
AC_LANG_POP([C++])

AC_CHECK_HEADERS([assert.h ctype.h fcntl.h inttypes.h math.h ncurses.h pthread.h semaphore.h stdio.h stdlib.h sys/mman.h sys/stat.h sys/types.h unistd.h ])

###
#AC_CHECK_HEADERS([boost/thread/thread.hpp],[],[AC_MSG_ERROR(You need the include Boost header to build $PACKAGE])],[])
###

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T



# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset regcomp sqrt strdup strncasecmp strtol])




AC_CONFIG_FILES([Makefile])
AC_OUTPUT
