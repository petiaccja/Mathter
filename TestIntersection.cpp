//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Geometry.hpp"


using namespace mathter;


TEST_CASE("Intersection - Line-plane 3D", "[Intersection]") {
	Hyperplane<float, 3> plane({ 1,2,3 }, Vector<float, 3>{ 3,2,1 }.Normalized());
	Line<float, 3> line({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized());

	auto intersection = Intersect(plane, line);
	REQUIRE(intersection.Intersecting());
	Vector<float, 3> point = intersection.Point();
	Vector<float, 3> expected(-14, 21, 10);
	float angle = 90 - acos(Vector<float, 3>::Dot(plane.Normal(), line.Direction()))*180.f / 3.1415926f;
	float distance = (point - line.Base()).Length();
	REQUIRE(point.Approx() == expected);
}


TEST_CASE("Intersection - LineSegment-plane 3D", "[Intersection]") {
	Hyperplane<float, 3> plane({ 1,2,3 }, Vector<float, 3>{ 3, 2, 1 }.Normalized());
	LineSegment<float, 3> succeedLine({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized(), 22.45f * 1.5f);
	LineSegment<float, 3> failLine({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized(), 22.45f * 0.8f);

	auto interSucceed = Intersect(plane, succeedLine);
	auto interFail = Intersect(plane, failLine);

	REQUIRE(interSucceed.Intersecting());
	float paramexp = 0.666666f;
	REQUIRE(interSucceed.InterpolParameter() == Approx(paramexp));
	REQUIRE(interSucceed.Point().Approx() == Vector<float, 3>(-14, 21, 10));

	REQUIRE_FALSE(interFail.Intersecting());
}


TEST_CASE("Intersection - Line-hyperplane", "[Intersection]") {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0 - 100,0 }, { 2 - 100,3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0 - 100,3 }, { 2 - 100,0 });

	line2 = Hyperplane<float, 2>(line1);

	Hyperplane<float, 2> plane1(line1);
	Hyperplane<float, 2> plane2(line2);

	REQUIRE(plane1.Normal().Approx() == plane2.Normal());
	REQUIRE(Approx(plane1.Scalar()) == plane2.Scalar());
}


TEST_CASE("Intersection - Line-line 2D", "[Intersection]") {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0-100,0 }, { 2-100,3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0-100,3 }, { 2-100,0 });

	auto inter = Intersect(line1, line2);
	auto point = inter.Point();

	REQUIRE(inter.Intersecting());
	REQUIRE(point.Approx() == Vector<float, 2>{ 1.0f-100, 1.5f });
}


TEST_CASE("Intersection - LineSegment-line", "[Intersection]") {
	Line<float, 2> line1 = LineSegment<float, 2>({ 0-100,0 }, { 2-100,3 }).Line();
	Line<float, 2> line2 = LineSegment<float, 2>({ 0-100,3 }, { 2-100,0 }).Line();

	auto inter = Intersect(line1, line2);

	REQUIRE(inter.Intersecting());
	REQUIRE(inter.Point().Approx() == Vector<float, 2>{ 1.0f-100, 1.5f });
}


TEST_CASE("Intersection - LineSegment-LineSegment", "[Intersection]") {
	LineSegment<float, 2> line1({ 0,0 }, { 2,3 });
	LineSegment<float, 2> lineSuc({ 0,4 }, { 2,-1 });
	LineSegment<float, 2> lineFail({ 0,1 }, { 2,4 });

	auto interSuc = Intersect(line1, lineSuc);
	auto interFail = Intersect(line1, lineFail);

	REQUIRE(interSuc.Intersecting());
	REQUIRE(interSuc.InterpolParameter1() == Approx(0.5f));
	REQUIRE(interSuc.InterpolParameter2() == Approx(0.5f));
	REQUIRE(interSuc.Point().Approx() == Vector<float, 2>{ 1.0f, 1.5f });

	REQUIRE_FALSE(interFail.Intersecting());
}