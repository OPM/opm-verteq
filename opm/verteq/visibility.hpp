#ifndef OPM_VERTEQ_VISIBILITY_HPP_INCLUDED
#define OPM_VERTEQ_VISIBILITY_HPP_INCLUDED

// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0

/* common visibility macros */
#ifndef OPM_VERTEQ_VISIBILITY_INCLUDED
#include <opm/verteq/utility/visibility.h>
#endif /* OPM_VERTEQ_VISIBILITY_INCLUDED */

/* special visibility macros for this module */
#if defined(opmverteq_EXPORTS)
#  define OPM_VERTEQ_PUBLIC  SYMBOL_IS_EXPORTED
#else
#  define OPM_VERTEQ_PUBLIC  SYMBOL_IS_IMPORTED
#endif
#define OPM_VERTEQ_PRIVATE  SYMBOL_IS_LOCALDEF

#endif /* OPM_VERTEQ_VISIBILITY_HPP_INCLUDED */
