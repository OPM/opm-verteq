# - Hide symbols in dynamic shared objects by default
#
# This give better encapsulation and link time.

include (AddOptions)
include (UseCompVer)
include (TestCXXAcceptsFlag)

get_gcc_version (CXX GCC_VERSION)
if (GCC_VERSION)
  if (GCC_VERSION VERSION_EQUAL 4.0 OR GCC_VERSION VERSION_GREATER 4.0)
	# unless we explicity mark classes and functions as exported, they
	# shouldn't be (like they are on Windows by default)
	check_cxx_accepts_flag ("-fvisibility=hidden" HAVE_VISIBILITY_HIDDEN)
	if (HAVE_VISIBILITY_HIDDEN)
	  add_options ("C;CXX" ALL_BUILDS "-fvisibility=hidden")
	endif (HAVE_VISIBILITY_HIDDEN)

	# mark inline members in the class as hidden (they will also be
	# inlined in the client code)
	check_cxx_accepts_flag ("-fvisibility-inlines-hidden" HAVE_INLINES_HIDDEN)
	if (HAVE_INLINES_HIDDEN)
	  add_options ("C;CXX" ALL_BUILDS "-fvisibility-inlines-hidden")
	endif (HAVE_INLINES_HIDDEN)
	
	# don't export symbols that come from static libraries
	check_cxx_accepts_flag ("-Wl,--exclude-libs,ALL" HAVE_EXCLUDE_LIBS)
	if (HAVE_EXCLUDE_LIBS)
	  add_options ("C;CXX" ALL_BUILDS "-Wl,--exclude-libs,ALL")
	endif (HAVE_EXCLUDE_LIBS)
  endif (GCC_VERSION VERSION_EQUAL 4.0 OR GCC_VERSION VERSION_GREATER 4.0)
endif (GCC_VERSION)
  
