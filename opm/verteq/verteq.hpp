#ifndef OPM_VERTEQ_INCLUDED
#define OPM_VERTEQ_INCLUDED

#include <opm/core/utility/visibility.h>

#if defined(opmverteq_EXPORTS)
#  define OPM_VERTEQ_PUBLIC  SYMBOL_IS_EXPORTED
#else
#  define OPM_VERTEQ_PUBLIC  SYMBOL_IS_IMPORTED
#endif
#define OPM_VERTEQ_PRIVATE  SYMBOL_IS_LOCALDEF


#endif /* OPM_VERTEQ_INCLUDED */
