dnl -*- autoconf -*-

dnl locate opm-verteq library itself; this macro is called by every module
dnl that depends on opm-verteq.
AC_DEFUN([OPM_VERTEQ_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-verteq],[1.0],[OPM Vertical Equilibrium Library])
])

dnl find all prerequisites of opm-verteq; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_VERTEQ_CHECKS],[])
