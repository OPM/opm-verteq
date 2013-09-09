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

#define BOOST_TEST_MODULE Navigation
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

// interface to module we are testing
#include <opm/verteq/nav.hpp>

BOOST_AUTO_TEST_SUITE (Navigation)

BOOST_AUTO_TEST_CASE (dir_equals)
{
	Dir dec = Dir::DEC;
	Dir inc = Dir::INC;
	Dir dec_twin = Dir::DEC;
	Dir inc_twin = Dir::INC;
	BOOST_REQUIRE_EQUAL (dec, dec);
	BOOST_REQUIRE_EQUAL (dec, dec_twin);
	BOOST_REQUIRE_EQUAL (inc, inc);
	BOOST_REQUIRE_EQUAL (inc, inc_twin);
	BOOST_REQUIRE_NE (dec, inc);
}

BOOST_AUTO_TEST_CASE (dim_equals)
{
	Dim2D x = Dim2D::X;
	Dim2D y = Dim2D::Y;
	Dim2D x_twin = Dim2D::X;
	Dim2D y_twin = Dim2D::Y;
	BOOST_REQUIRE_EQUAL (x, x);
	BOOST_REQUIRE_EQUAL (x, x_twin);
	BOOST_REQUIRE_EQUAL (y, y);
	BOOST_REQUIRE_EQUAL (y, y_twin);
	BOOST_REQUIRE_NE (x, y);
}

BOOST_AUTO_TEST_CASE (dir_opposite)
{
	Dir dec = Dir::DEC;
	Dir inc = Dir::INC;
	BOOST_REQUIRE_EQUAL (dec.opposite (), inc);
	BOOST_REQUIRE_EQUAL (inc.opposite (), dec);
}

BOOST_AUTO_TEST_CASE (dim_orthogonal)
{
	Dim2D x = Dim2D::X;
	Dim2D y = Dim2D::Y;
	BOOST_REQUIRE_EQUAL (x.orthogonal (), y);
	BOOST_REQUIRE_EQUAL (y.orthogonal (), x);
}

BOOST_AUTO_TEST_CASE (facetag)
{
	// see opm/core/grid.h, UnstructuredGrid::face_tag
	BOOST_REQUIRE_EQUAL (Side3D (Dim2D::X, Dir::DEC).facetag (), 0);
	BOOST_REQUIRE_EQUAL (Side3D (Dim2D::X, Dir::INC).facetag (), 1);
	BOOST_REQUIRE_EQUAL (Side3D (Dim2D::Y, Dir::DEC).facetag (), 2);
	BOOST_REQUIRE_EQUAL (Side3D (Dim2D::Y, Dir::INC).facetag (), 3);
	BOOST_REQUIRE_EQUAL (Side3D (Dim3D::Z, Dir::DEC).facetag (), 4);
	BOOST_REQUIRE_EQUAL (Side3D (Dim3D::Z, Dir::INC).facetag (), 5);
}

BOOST_AUTO_TEST_CASE (fromtag)
{
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (0), Side3D (Dim2D::X, Dir::DEC));
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (1), Side3D (Dim2D::X, Dir::INC));
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (2), Side3D (Dim2D::Y, Dir::DEC));
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (3), Side3D (Dim2D::Y, Dir::INC));
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (4), Side3D (Dim3D::Z, Dir::DEC));
	BOOST_REQUIRE_EQUAL (Side3D::from_tag (5), Side3D (Dim3D::Z, Dir::INC));
}

BOOST_AUTO_TEST_CASE (face_enum)
{
	std::vector <int> num (Side3D::COUNT, 0);
	std::vector <int> ones (Side3D::COUNT, 1);

	// make sure that each side it hit once and only once
	for (const Side3D* s = Side3D::begin (); s != Side3D::end (); ++s) {
		num[s->facetag()]++;
	}
	BOOST_CHECK_EQUAL_COLLECTIONS (num.begin(), num.end(),
	                               ones.begin(), ones.end());
}

BOOST_AUTO_TEST_CASE (pivot)
{
	Corn3D lft (Dir::DEC, Dir::DEC, Dir::DEC); // left, front, top
	Corn3D rft (Dir::INC, Dir::DEC, Dir::DEC); // right, front, top
	Corn3D lbt (Dir::DEC, Dir::INC, Dir::DEC); // left, back, top
	Corn3D rbt (Dir::INC, Dir::INC, Dir::DEC); // right, back, top
	Corn3D lfd (Dir::DEC, Dir::DEC, Dir::INC); // left, front, down
	Corn3D rfd (Dir::INC, Dir::DEC, Dir::INC); // right, front, down
	Corn3D lbd (Dir::DEC, Dir::INC, Dir::INC); // left, back, down
	Corn3D rbd (Dir::INC, Dir::INC, Dir::INC); // right, back, down

	BOOST_REQUIRE_EQUAL (lft.pivot (Dim2D::X, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (lft.pivot (Dim2D::X, Dir::INC), rft);
	BOOST_REQUIRE_EQUAL (lft.pivot (Dim2D::Y, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (lft.pivot (Dim2D::Y, Dir::INC), lbt);
	BOOST_REQUIRE_EQUAL (lft.pivot (Dim3D::Z, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (lft.pivot (Dim3D::Z, Dir::INC), lfd);

	BOOST_REQUIRE_EQUAL (rft.pivot (Dim2D::X, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (rft.pivot (Dim2D::X, Dir::INC), rft);
	BOOST_REQUIRE_EQUAL (rft.pivot (Dim2D::Y, Dir::DEC), rft);
	BOOST_REQUIRE_EQUAL (rft.pivot (Dim2D::Y, Dir::INC), rbt);
	BOOST_REQUIRE_EQUAL (rft.pivot (Dim3D::Z, Dir::DEC), rft);
	BOOST_REQUIRE_EQUAL (rft.pivot (Dim3D::Z, Dir::INC), rfd);

	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim2D::X, Dir::DEC), lbt);
	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim2D::X, Dir::INC), rbt);
	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim2D::Y, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim2D::Y, Dir::INC), lbt);
	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim3D::Z, Dir::DEC), lbt);
	BOOST_REQUIRE_EQUAL (lbt.pivot (Dim3D::Z, Dir::INC), lbd);

	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim2D::X, Dir::DEC), lbt);
	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim2D::X, Dir::INC), rbt);
	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim2D::Y, Dir::DEC), rft);
	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim2D::Y, Dir::INC), rbt);
	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim3D::Z, Dir::DEC), rbt);
	BOOST_REQUIRE_EQUAL (rbt.pivot (Dim3D::Z, Dir::INC), rbd);

	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim2D::X, Dir::DEC), lfd);
	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim2D::X, Dir::INC), rfd);
	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim2D::Y, Dir::DEC), lfd);
	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim2D::Y, Dir::INC), lbd);
	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim3D::Z, Dir::DEC), lft);
	BOOST_REQUIRE_EQUAL (lfd.pivot (Dim3D::Z, Dir::INC), lfd);

	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim2D::X, Dir::DEC), lfd);
	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim2D::X, Dir::INC), rfd);
	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim2D::Y, Dir::DEC), rfd);
	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim2D::Y, Dir::INC), rbd);
	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim3D::Z, Dir::DEC), rft);
	BOOST_REQUIRE_EQUAL (rfd.pivot (Dim3D::Z, Dir::INC), rfd);

	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim2D::X, Dir::DEC), lbd);
	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim2D::X, Dir::INC), rbd);
	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim2D::Y, Dir::DEC), lfd);
	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim2D::Y, Dir::INC), lbd);
	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim3D::Z, Dir::DEC), lbt);
	BOOST_REQUIRE_EQUAL (lbd.pivot (Dim3D::Z, Dir::INC), lbd);

	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim2D::X, Dir::DEC), lbd);
	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim2D::X, Dir::INC), rbd);
	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim2D::Y, Dir::DEC), rfd);
	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim2D::Y, Dir::INC), rbd);
	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim3D::Z, Dir::DEC), rbt);
	BOOST_REQUIRE_EQUAL (rbd.pivot (Dim3D::Z, Dir::INC), rbd);
}

BOOST_AUTO_TEST_CASE (cart_elems)
{
	Cart2D c (2, 2);
	Coord2D _00 (0, 0); Coord2D _10 (1, 0);
	Coord2D _01 (0, 1); Coord2D _11 (1, 1);
	BOOST_REQUIRE_EQUAL (c.num_elems (), 4);
	BOOST_REQUIRE_EQUAL (c.cart_ndx (_00), 0);
	BOOST_REQUIRE_EQUAL (c.cart_ndx (_10), 1);
	BOOST_REQUIRE_EQUAL (c.cart_ndx (_01), 2);
	BOOST_REQUIRE_EQUAL (c.cart_ndx (_11), 3);

	BOOST_REQUIRE_EQUAL (c.coord (0), _00);
	BOOST_REQUIRE_EQUAL (c.coord (1), _10);
	BOOST_REQUIRE_EQUAL (c.coord (2), _01);
	BOOST_REQUIRE_EQUAL (c.coord (3), _11);
}

BOOST_AUTO_TEST_CASE (cart_nodes)
{
	Cart2D c (2, 2);
	Coord2D _00 (0, 0); Coord2D _10 (1, 0);
	Coord2D _01 (0, 1); Coord2D _11 (1, 1);
	Corn2D lf (Dir::DEC, Dir::DEC); // left, front
	Corn2D rf (Dir::INC, Dir::DEC); // right, front
	Corn2D lb (Dir::DEC, Dir::INC); // left, back
	Corn2D rb (Dir::INC, Dir::INC); // right, back

	BOOST_REQUIRE_EQUAL (c.num_nodes (), 9);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_00, lf), 0);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_00, rf), 1);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_00, lb), 3);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_00, rb), 4);

	BOOST_REQUIRE_EQUAL (c.node_ndx (_10, lf), 1);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_10, rf), 2);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_10, lb), 4);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_10, rb), 5);

	BOOST_REQUIRE_EQUAL (c.node_ndx (_01, lf), 3);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_01, rf), 4);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_01, lb), 6);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_01, rb), 7);

	BOOST_REQUIRE_EQUAL (c.node_ndx (_11, lf), 4);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_11, rf), 5);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_11, lb), 7);
	BOOST_REQUIRE_EQUAL (c.node_ndx (_11, rb), 8);
}

BOOST_AUTO_TEST_CASE (cart_faces)
{
	Cart2D c (2, 2);
	Coord2D _00 (0, 0); Coord2D _10 (1, 0);
	Coord2D _01 (0, 1); Coord2D _11 (1, 1);

	// notice the terminology; the "left" side is the one
	// which is in decreasing X direction from the center
	Side2D l (Dim2D::X, Dir::DEC); // left
	Side2D r (Dim2D::X, Dir::INC); // right
	Side2D f (Dim2D::Y, Dir::DEC); // front
	Side2D b (Dim2D::Y, Dir::INC); // back

	BOOST_REQUIRE_EQUAL (c.num_faces (), 12);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_00, l),  2);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_00, r),  3);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_00, f),  0);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_00, b),  5);

	BOOST_REQUIRE_EQUAL (c.face_ndx (_10, l),  3);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_10, r),  4);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_10, f),  1);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_10, b),  6);

	BOOST_REQUIRE_EQUAL (c.face_ndx (_01, l),  7);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_01, r),  8);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_01, f),  5);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_01, b), 10);

	BOOST_REQUIRE_EQUAL (c.face_ndx (_11, l),  8);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_11, r),  9);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_11, f),  6);
	BOOST_REQUIRE_EQUAL (c.face_ndx (_11, b), 11);
}

BOOST_AUTO_TEST_SUITE_END ()
