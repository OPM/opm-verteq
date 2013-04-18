// Copyright (C) 2013 Uni Research AS
// This file is licensed under the GNU General Public License v3.0
#include <opm/verteq/nav.hpp>

const Dir Dir::DEC (0);
const Dir Dir::INC (1);

const Dim2D Dim2D::X (0);
const Dim2D Dim2D::Y (1);
const Dim3D Dim3D::Z (2);

template <typename Dim> Side <Dim>
Side <Dim>::from_tag (int tag) {
	// direction is the minor bits in the enumeration
	const div_t bits = div (tag, Dir::COUNT);
	const int dir_val = bits.rem;
	const int dim_val = bits.quot;
	return Side <Dim> (Dim (dim_val), Dir (dir_val));
}

// template instantiation to satisfy linker
template Side <Dim2D> Side <Dim2D>::from_tag (int);
template Side <Dim3D> Side <Dim3D>::from_tag (int);

// these needs to be initialized here instead of in the header
// because vector::resize takes a reference to the data and not
// a value as a parameter (to avoid copying)
const int Cart2D::NO_ELEM = -1;
const int Cart2D::NO_FACE = -1;
const int Cart2D::NO_NODE = -1;
