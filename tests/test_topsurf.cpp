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

#define BOOST_TEST_MODULE TopSurfTest
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

// interface to module we are testing
#include <opm/verteq/topsurf.hpp>

// utility modules (to setup fine grid)
#include <opm/core/grid/cart_grid.h>

/**
 * Sample grid used to manually test the top-surface generation.
 *
 * The numbers in this test are intentionally hardcoded. If we change
 * the system of generation, for instance the orientation of the faces
 * around each cell, then these tests will fail and will have to be
 * updated if the change was planned.
 *
 * The top surface looks like this (using V prefix for vertices, F
 * for faces, and E for elements):
 *
 * 1.0 -    V4------F3-----V5-----F6-----V6-----F9-----V7
 *     |     |             |             |             |
 *     |     |             |             |             |
 * 0.5 -     F0     E0     F1     E1     F4      E3    F7
 *     |     |             |             |             |
 *     |     |             |             |             |
 * 0.0 -    V0-----F2------V1-----F5-----V2-----F8-----V3
 *     |
 *     |-----|------|------|------|------|------|------|
 *          0.0    0.5    1.0    1.5    2.0    2.5    3.0
 */
struct HardcodedGrids {
	UnstructuredGrid* g; // fine grid
	Opm::TopSurf* ts;    // coarse grid

	HardcodedGrids () {
		g = create_grid_cart3d (3, 1, 2);
		ts = Opm::TopSurf::create (*g);
	}

	~HardcodedGrids () {
		// ts has a destructor!
		delete ts;
		destroy_grid (g);
	}
};

BOOST_FIXTURE_TEST_SUITE (TopSurfTest, HardcodedGrids)

BOOST_AUTO_TEST_CASE (format)
{
	// original grid was 3D, top surface should be 2D
	BOOST_REQUIRE_EQUAL (ts->dimensions, 2);
	BOOST_REQUIRE_EQUAL (ts->number_of_cells,  3);
	BOOST_REQUIRE_EQUAL (ts->number_of_faces, 10);
	BOOST_REQUIRE_EQUAL (ts->number_of_nodes,  8);

	// top-surface should be Cartesian, too
	const int cartdims[] = { 3, 1 };
	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cartdims,
	                               ts->cartdims+ts->dimensions,
	                               cartdims,
	                               cartdims+sizeof(cartdims)/sizeof(cartdims[0]));

	const int global_cell[] = { 0, 1, 2 };
	BOOST_CHECK_EQUAL_COLLECTIONS (ts->global_cell,
	                               ts->global_cell+ts->number_of_cells,
	                               global_cell,
	                               global_cell+sizeof(global_cell)/sizeof(global_cell[0]));
}

BOOST_AUTO_TEST_CASE (embedding)
{
	const double vertex_coords[] = {
	  0., 0.,    1., 0.,    2., 0.,    3., 0.,
	  0., 1.,    1., 1.,    2., 1.,    3., 1.,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->node_coordinates,
	                               ts->node_coordinates+ts->dimensions*ts->number_of_nodes,
	                               vertex_coords,
	                               vertex_coords+sizeof(vertex_coords)/sizeof(vertex_coords[0]));

	const double element_coords[] = {
	  0.5, 0.5,    1.5, 0.5,    2.5, 0.5,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cell_centroids,
	                               ts->cell_centroids+ts->dimensions*ts->number_of_cells,
	                               element_coords,
	                               element_coords+sizeof(element_coords)/sizeof(element_coords[0]));

	const double face_coords[] = {
	  0.0, 0.5,
	  1.0, 0.5,    0.5, 0.0,    0.5, 1.0,
	  2.0, 0.5,    1.5, 0.0,    1.5, 1.0,
	  3.0, 0.5,    2.5, 0.0,    2.5, 1.0,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_centroids,
	                               ts->face_centroids+ts->dimensions*ts->number_of_faces,
	                               face_coords,
	                               face_coords+sizeof(face_coords)/sizeof(face_coords[0]));
}

BOOST_AUTO_TEST_CASE (orientation)
{
	const double face_normals[] = {
	  1.0, 0.0,
	  1.0, 0.0,    0.0, 1.0,    0.0, 1.0,
	  1.0, 0.0,    0.0, 1.0,    0.0, 1.0,
	  1.0, 0.0,    0.0, 1.0,    0.0, 1.0,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_normals,
	                               ts->face_normals+ts->dimensions*ts->number_of_faces,
	                               face_normals,
	                               face_normals+sizeof(face_normals)/sizeof(face_normals[0]));

	const int facetag[] = {
	  0,    1,    2,    3,
	  0,    1,    2,    3,
	  0,    1,    2,    3,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cell_facetag,
	                               ts->cell_facetag+ts->cell_facepos[ts->number_of_cells],
	                               facetag,
	                               facetag+sizeof(facetag)/sizeof(facetag[0]));
}

BOOST_AUTO_TEST_CASE (sizes)
{
	const double areas[] = {
	  1.,
	  1.,    1.,    1.,
	  1.,    1.,    1.,
	  1.,    1.,    1.,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_areas,
	                               ts->face_areas+ts->number_of_faces,
	                               areas,
	                               areas+sizeof(areas)/sizeof(areas[0]));

	const double volumes[] = {
	  1.,    1.,    1.,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cell_volumes,
	                               ts->cell_volumes+ts->number_of_cells,
	                               volumes,
	                               volumes+sizeof(volumes)/sizeof(volumes[0]));
}

BOOST_AUTO_TEST_CASE (topology)
{
	// each edge in 2D always have two and only two endpoints
	const int nodepos[] = {
	  0,
	  2,    4,    6,
	  8,   10,   12,
	 14,   16,   18,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_nodepos,
	                               ts->face_nodepos+ts->number_of_faces,
	                               nodepos,
	                               nodepos+sizeof(nodepos)/sizeof(nodepos[0]));

	const int nodes[] = {
	  0, 4,
	  1, 5,    1, 0,    5, 4,
	  2, 6,    2, 1,    6, 5,
	  3, 7,    3, 2,    7, 6,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_nodes,
	                               ts->face_nodes+ts->face_nodepos[ts->number_of_faces],
	                               nodes,
	                               nodes+sizeof(nodes)/sizeof(nodes[0]));

	// all elements are quadrilaterals
	const int facepos[] = {
	  0,    4,    8,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cell_facepos,
	                               ts->cell_facepos+ts->number_of_cells,
	                               facepos,
	                               facepos+sizeof(facepos)/sizeof(facepos[0]));

	const int faces[] = {
	  0,    1,     2,     3,
	  1,    4,     5,     6,
	  4,    7,     8,     9,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->cell_faces,
	                               ts->cell_faces+ts->cell_facepos[ts->number_of_cells],
	                               faces,
	                               faces+sizeof(faces)/sizeof(faces[0]));
}

BOOST_AUTO_TEST_CASE (neighbours)
{
	const int face_cells[] = {
	  -1,  0,
	   0,  1,    -1, 0,    0, -1,
	   1,  2,    -1, 1,    1, -1,
	   2, -1,    -1, 2,    2, -1,
	};

	BOOST_CHECK_EQUAL_COLLECTIONS (ts->face_cells,
	                               ts->face_cells+2*ts->number_of_faces,
	                               face_cells,
	                               face_cells+sizeof(face_cells)/sizeof(face_cells[0]));
}

BOOST_AUTO_TEST_SUITE_END ()
