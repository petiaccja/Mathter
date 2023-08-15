// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Geometry.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>


using namespace mathter;
using Catch::Approx;


using Ray3 = Ray<float, 3>;
using Vec3 = Vector<float, 3>;
using Vec4 = Vector<float, 4>;
using Triangle = Triangle3D<float>;


TEST_CASE("Intersection - Line-plane 3D", "[Intersection]") {
	Hyperplane<float, 3> plane({ 1, 2, 3 }, Normalize(Vector<float, 3>{ 3, 2, 1 }));
	Line<float, 3> line({ -2, 3, 4 }, Normalize(Vector<float, 3>{ -2, 3, 1 }));

	auto intersection = Intersect(plane, line);
	REQUIRE(intersection.Intersecting());
	Vector<float, 3> point = intersection.Point();
	Vector<float, 3> expected(-14, 21, 10);
	float angle = 90 - acos(Dot(plane.Normal(), line.Direction())) * 180.f / 3.1415926f;
	float distance = Length(point - line.Base());
	REQUIRE(ApproxVec(point) == expected);
}


TEST_CASE("Intersection - LineSegment-plane 3D", "[Intersection]") {
	Hyperplane<float, 3> plane({ 1, 2, 3 }, Normalize(Vector<float, 3>{ 3, 2, 1 }));
	LineSegment<float, 3> succeedLine({ -2, 3, 4 }, Normalize(Vector<float, 3>{ -2, 3, 1 }), 22.45f * 1.5f);
	LineSegment<float, 3> failLine({ -2, 3, 4 }, Normalize(Vector<float, 3>{ -2, 3, 1 }), 22.45f * 0.8f);

	auto interSucceed = Intersect(plane, succeedLine);
	auto interFail = Intersect(plane, failLine);

	REQUIRE(interSucceed.Intersecting());
	float paramexp = 0.666666f;
	REQUIRE(interSucceed.InterpolParameter() == Approx(paramexp));
	REQUIRE(ApproxVec(interSucceed.Point()) == Vector<float, 3>(-14, 21, 10));

	REQUIRE_FALSE(interFail.Intersecting());
}


TEST_CASE("Intersection - Line-hyperplane", "[Intersection]") {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0 - 100, 0 }, { 2 - 100, 3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0 - 100, 3 }, { 2 - 100, 0 });

	line2 = Hyperplane<float, 2>(line1);

	Hyperplane<float, 2> plane1(line1);
	Hyperplane<float, 2> plane2(line2);

	REQUIRE(ApproxVec(plane1.Normal()) == plane2.Normal());
	REQUIRE(Approx(plane1.Scalar()) == plane2.Scalar());
}


TEST_CASE("Intersection - Line-line 2D", "[Intersection]") {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0 - 100, 0 }, { 2 - 100, 3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0 - 100, 3 }, { 2 - 100, 0 });

	auto inter = Intersect(line1, line2);
	auto point = inter.Point();

	REQUIRE(inter.Intersecting());
	REQUIRE(ApproxVec(point) == Vector<float, 2>{ 1.0f - 100, 1.5f });
}


TEST_CASE("Intersection - LineSegment-line", "[Intersection]") {
	Line<float, 2> line1 = LineSegment<float, 2>({ 0 - 100, 0 }, { 2 - 100, 3 }).Line();
	Line<float, 2> line2 = LineSegment<float, 2>({ 0 - 100, 3 }, { 2 - 100, 0 }).Line();

	auto inter = Intersect(line1, line2);

	REQUIRE(inter.Intersecting());
	REQUIRE(ApproxVec(inter.Point()) == Vector<float, 2>{ 1.0f - 100, 1.5f });
}


TEST_CASE("Intersection - LineSegment-LineSegment", "[Intersection]") {
	LineSegment<float, 2> line1({ 0, 0 }, { 2, 3 });
	LineSegment<float, 2> lineSuc({ 0, 4 }, { 2, -1 });
	LineSegment<float, 2> lineFail1({ 0, 1 }, { 2, 4 });
	LineSegment<float, 2> lineFail2({ 0, 2 }, { 3, 4 });

	auto interSuc = Intersect(line1, lineSuc);
	auto interFail1 = Intersect(line1, lineFail1);
	auto interFail2 = Intersect(line1, lineFail2);

	REQUIRE(interSuc.Intersecting());
	REQUIRE(interSuc.InterpolParameter1() == Approx(0.5f));
	REQUIRE(interSuc.InterpolParameter2() == Approx(0.5f));
	REQUIRE(ApproxVec(interSuc.Point()) == Vector<float, 2>{ 1.0f, 1.5f });

	REQUIRE_FALSE(interFail1.Intersecting());
	REQUIRE_FALSE(interFail2.Intersecting());
}


TEST_CASE("Ray - construction", "[Ray-triangle intersect]") {
	Ray3 ray{ Vec3(1, 2, 3), Normalize(Vec3(2, 4, 6)) };
	REQUIRE(ApproxVec(ray.Direction()) == Normalize(Vec3(2, 4, 6)));
	REQUIRE(ApproxVec(ray.Base()) == Vec3(1, 2, 3));
}


TEST_CASE("Ray-tri - aligned hit", "[Ray-triangle intersect]") {
	Ray3 ray{ Vec3(0.5f, 0.0f, 0.5f), Vec3(0, 1, 0) };
	Triangle tri{ { 0, 1, 0 }, { 0.5f, 1, 1 }, { 1, 1, 0 } };

	auto intersection = Intersect(ray, tri);
	REQUIRE(intersection.IsIntersecting() == true);
	REQUIRE(ApproxVec(intersection.Point()) == Vec3(0.5f, 1.0f, 0.5f));
}


TEST_CASE("Ray-tri - aligned miss", "[Ray-triangle intersect]") {
	Ray3 ray{ Vec3(1.5f, 0.0f, 0.5f), Vec3(0, 1, 0) };
	Triangle tri{ { 0, 1, 0 }, { 0.5f, 1, 1 }, { 1, 1, 0 } };

	auto intersection = Intersect(ray, tri);
	REQUIRE(intersection.IsIntersecting() == false);
}


TEST_CASE("Ray-tri - interpol position", "[Ray-triangle intersect]") {
	Ray3 ray{ Vec3(0.35f, 0.0f, 0.55f), Vec3(0, 1, 0) };
	Triangle tri{ { 0, 1, 0 }, { 0.5f, 1, 1 }, { 1, 1, 0 } };

	auto intersection = Intersect(ray, tri);
	REQUIRE(intersection.IsIntersecting() == true);

	Vec3 pxpos = intersection.Interpolate(tri.a, tri.b, tri.c);
	REQUIRE(ApproxVec(pxpos) == intersection.Point());
}


TEST_CASE("Ray-tri - interpol property", "[Ray-triangle intersect]") {
	Ray3 ray{ Vec3(0.45f, 0.0f, 0.55f), Vec3(0, 1, 0) };
	Triangle tri{ { 0, 1, 0 }, { 0.5f, 1, 1 }, { 1, 1, 0 } };

	auto intersection = Intersect(ray, tri);
	REQUIRE(intersection.IsIntersecting() == true);

	Vec4 colors[3] = {
		{ 1, 0, 0, 1 },
		{ 0, 1, 0, 1 },
		{ 0, 0, 1, 1 },
	};

	Vec4 pxcolor = intersection.Interpolate(colors[0], colors[1], colors[2]);
	REQUIRE(pxcolor.y > pxcolor.z);
	REQUIRE(pxcolor.x > pxcolor.z);
	REQUIRE(pxcolor.y > pxcolor.x);
	REQUIRE(Approx(pxcolor.w) == 1);
}
