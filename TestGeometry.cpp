//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Geometry.hpp"

using namespace mathter;


TEST_CASE("Bezier curve", "[Geometry]") {
	Vector<float, 2> p1, p2, p3;
	p1 = {0,0};
	p2 = {1,2};
	p3 = {2,0};
    BezierCurve<float, 2, 2> curve = { p1, p2, p3 };

	auto point = curve(0.5f);
	Vector<float, 2> exp = {1, 1};
	REQUIRE(point.Approx() == exp);
}