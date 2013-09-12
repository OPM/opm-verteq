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

#define INVALID -1

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
	int pos[4];
	int data[8];

	Opm::rlw_int m;

	EncodedSparseMatrix () :
	  num (2)
	, pos { 0,
	        3,
	        7,
	        INVALID }
	, data { 11, 12, 13,
	         21, 22, 23, 24,
	         INVALID }
	, m (num, pos, data)
	{}
};

BOOST_FIXTURE_TEST_SUITE (RunLen, EncodedSparseMatrix)

BOOST_AUTO_TEST_CASE (operator_index)
{
	// column indexing
	BOOST_REQUIRE_EQUAL (m[0], &data[0]);
	BOOST_REQUIRE_EQUAL (m[1], &data[3]);
}

BOOST_AUTO_TEST_CASE (cols)
{
	// number of columns
	BOOST_REQUIRE_EQUAL (m.cols(), 2);
}

BOOST_AUTO_TEST_CASE (size)
{
	// size of each column
	BOOST_REQUIRE_EQUAL (m.size (0), 3);
	BOOST_REQUIRE_EQUAL (m.size (1), 4);
}

BOOST_AUTO_TEST_CASE (first)
{
	// first element in each column
	BOOST_REQUIRE_EQUAL (m[0][0], 11);
	BOOST_REQUIRE_EQUAL (m[1][0], 21);
}

BOOST_AUTO_TEST_CASE (last)
{
	// last item in each column
	BOOST_REQUIRE_EQUAL (m.last (0), 13);
	BOOST_REQUIRE_EQUAL (m.last (1), 24);
}

BOOST_AUTO_TEST_CASE (copy)
{
	// copy constructor
	Opm::rlw_int a (m);
	BOOST_REQUIRE_EQUAL (m.cols(), a.cols());
	for (int i = 0; i < m.cols(); ++i) {
		BOOST_REQUIRE_EQUAL (m[0], a[0]);
		BOOST_REQUIRE_EQUAL (m.size (i), a.size (i));
		for (int j = 0; j < m.size (i); ++j) {
			BOOST_REQUIRE_EQUAL (m[i][j], a[i][j]);
		}
	}
}

BOOST_AUTO_TEST_CASE (deep)
{
	// deep copy
	Opm::RunLenData <int> a (num, pos);
	BOOST_REQUIRE_EQUAL (m.cols (), a.cols ());
	for (int i = 0; i < m.cols(); ++i) {
		BOOST_REQUIRE_EQUAL (m.size (i), a.size (i));
		for (int j = 0; j < m.size (i); ++j) {
			a[i][j] = m[i][j];
		}
	}
	// read back
	for (int i = 0; i < a.cols(); ++i) {
		for (int j = 0; j < a.size (i); ++j) {
			BOOST_REQUIRE_EQUAL (a[i][j], m[i][j]);
			BOOST_REQUIRE_NE (a[i][j], INVALID);
		}
	}
}

BOOST_AUTO_TEST_SUITE_END ()
