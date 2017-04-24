#pragma warning(disable: 4244)
#include <gtest\gtest.h>

#include "Mathter/Geometry.hpp"


using namespace mathter;


TEST(Intersect, Line_Plane_3D) {
	Hyperplane<float, 3> plane({ 1,2,3 }, Vector<float, 3>{ 3,2,1 }.Normalized());
	Line<float, 3> line({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized());

	auto intersection = Intersect(plane, line);
	ASSERT_TRUE(intersection.Intersecting());
	Vector<float, 3> point = intersection.Point();
	Vector<float, 3> expected(-14, 21, 10);
	float angle = 90 - acos(Vector<float, 3>::Dot(plane.Normal(), line.Direction()))*180.f / 3.1415926f;
	float distance = (point - line.Base()).Length();
	ASSERT_TRUE(point.AlmostEqual(expected));
}


TEST(Intersect, LineSegment_Plane_3D) {
	Hyperplane<float, 3> plane({ 1,2,3 }, Vector<float, 3>{ 3, 2, 1 }.Normalized());
	LineSegment<float, 3> succeedLine({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized(), 22.45f * 1.5f);
	LineSegment<float, 3> failLine({ -2, 3, 4 }, Vector<float, 3>{ -2, 3, 1 }.Normalized(), 22.45f * 0.8f);

	auto interSucceed = Intersect(plane, succeedLine);
	auto interFail = Intersect(plane, failLine);

	ASSERT_TRUE(interSucceed.Intersecting());
	float paramexp = 0.666666f;
	ASSERT_NEAR(interSucceed.InterpolParameter(), paramexp, 0.001f);
	ASSERT_TRUE(interSucceed.Point().AlmostEqual(Vector<float, 3>(-14, 21, 10)));

	ASSERT_FALSE(interFail.Intersecting());
}


TEST(Intersect, Line_Hyperplane) {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0 - 100,0 }, { 2 - 100,3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0 - 100,3 }, { 2 - 100,0 });

	line2 = Hyperplane<float, 2>(line1);

	Hyperplane<float, 2> plane1(line1);
	Hyperplane<float, 2> plane2(line2);

	ASSERT_TRUE(plane1.Normal().AlmostEqual(plane2.Normal()));
	ASSERT_FLOAT_EQ(plane1.Scalar(), plane2.Scalar());
}


TEST(Intersect, Line_Line_2D) {
	Line<float, 2> line1 = Line<float, 2>::Through({ 0-100,0 }, { 2-100,3 });
	Line<float, 2> line2 = Line<float, 2>::Through({ 0-100,3 }, { 2-100,0 });

	auto inter = Intersect(line1, line2);
	auto point = inter.Point();

	ASSERT_TRUE(inter.Intersecting());
	ASSERT_TRUE(point.AlmostEqual({ 1.0f-100, 1.5f }));
}


TEST(Intersect, LineSegment_To_Line) {
	Line<float, 2> line1 = LineSegment<float, 2>({ 0-100,0 }, { 2-100,3 }).Line();
	Line<float, 2> line2 = LineSegment<float, 2>({ 0-100,3 }, { 2-100,0 }).Line();

	auto inter = Intersect(line1, line2);

	ASSERT_TRUE(inter.Intersecting());
	ASSERT_TRUE(inter.Point().AlmostEqual({ 1.0f-100, 1.5f }));
}


TEST(Intersect, LineSegment_LineSegment_2D) {
	LineSegment<float, 2> line1({ 0,0 }, { 2,3 });
	LineSegment<float, 2> lineSuc({ 0,4 }, { 2,-1 });
	LineSegment<float, 2> lineFail({ 0,1 }, { 2,4 });

	auto interSuc = Intersect(line1, lineSuc);
	auto interFail = Intersect(line1, lineFail);

	ASSERT_TRUE(interSuc.Intersecting());
	ASSERT_FLOAT_EQ(interSuc.InterpolParameter1(), 0.5f);
	ASSERT_FLOAT_EQ(interSuc.InterpolParameter2(), 0.5f);
	ASSERT_TRUE(interSuc.Point().AlmostEqual({ 1.0f, 1.5f }));

	ASSERT_FALSE(interFail.Intersecting());
}