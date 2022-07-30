//L=============================================================================
//L This software is distributed under the MIT license.
//L Copyright 2021 Péter Kardos
//L=============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include <Mathter/Geometry.hpp>
#include <Mathter/Common/Approx.hpp>

using namespace mathter;


TEST_CASE("Bezier curve", "[Geometry]") {
	Vector<float, 2> p1, p2, p3;
	p1 = {0,0};
	p2 = {1,2};
	p3 = {2,0};
    BezierCurve<float, 2, 2> curve = { p1, p2, p3 };

	auto point = curve(0.5f);
	Vector<float, 2> exp = {1, 1};
	REQUIRE(ApproxVec(point) == exp);
}