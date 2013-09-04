#ifdef HAVE_CONFIG_H
#  if HAVE_CONFIG_H
#    include <config.h>
#  endif
#endif /* HAVE_CONFIG_H */
#ifdef HAVE_DYNAMIC_BOOST_TEST
#  if HAVE_DYNAMIC_BOOST_TEST
#    define BOOST_TEST_DYN_LINK
#  endif
#endif /* HAVE_DYNAMIC_BOOST_TEST */

#define BOOST_TEST_MODULE RunLen
#include <boost/test/unit_test.hpp>

// interface to module we are testing
#include <opm/verteq/utility/runlen.hpp>

/**
 * Encoded sparse matrix of test data. It models this matrix:
 *
 * @code{.txt}
 *     ( 11 12 13 -   )
 *     ( 21 22 23 24  )
 * @endcode
 */
struct EncodedSparseMatrix {
	int num;
	int pos[3];
	int data[7];

	Opm::rlw_int m;

	EncodedSparseMatrix () :
	  num (2)
	, pos { 0,
	        3,
	        7 }
	, data { 11, 12, 13,
	         21, 22, 23, 24 }
	, m (num, pos, data)
	{}
};

BOOST_FIXTURE_TEST_SUITE (RunLen, EncodedSparseMatrix)

BOOST_AUTO_TEST_CASE(foo)
{
}

BOOST_AUTO_TEST_SUITE_END()
