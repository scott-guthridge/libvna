#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
AC_PREREQ([2.69])
AC_INIT([libvna], [0.1.0], [bugs@rompromity.net])
AC_CONFIG_SRCDIR([src/vnaconv-example.c])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CC_STDC
AC_PROG_INSTALL
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([yaml], [yaml_document_initialize])

# Checks for header files.
AC_CHECK_HEADERS([search.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([insque isascii mkdir random remque strcasecmp strdup vasprintf])

# Init libtool
AM_PROG_AR
LT_INIT

# Generate files.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
    Makefile
    rpm/Makefile
    rpm/libvna.spec
    src/Makefile
    src/tests/Makefile
])
AC_OUTPUT
