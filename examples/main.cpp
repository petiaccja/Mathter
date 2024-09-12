#include "Camera.hpp"
#include "Object.hpp"

#include "Mathter/Transforms.hpp"
#include "Mathter/Utility.hpp"

#include <array>
#include <iostream>
#include <thread>

using namespace mathter;

// The "texture" into which we render our image.
// We render to STDOUT, with characters instead of pixels.
using BackBuffer = Matrix<char, 25, 50>;


/// <summary> Apply a transform to all corners of a triangle. </summary>
Triangle3f TransformTriangle(const Triangle3f& triangle, const Mat44f& transform) {
	return Triangle3f(triangle.corners[0] * transform,
					  triangle.corners[1] * transform,
					  triangle.corners[2] * transform);
}


/// <summary> Determine if the triangle is front-facing, based on its winding order. </summary>
bool IsFrontFacing(const Triangle3f& triangle) {
	const auto edge1 = triangle.corners[1] - triangle.corners[0];
	const auto edge2 = triangle.corners[2] - triangle.corners[0];
	const auto normal = Cross(edge1, edge2);
	return normal.z <= 0.0f;
}


/// <summary> Rasterize an edge between <paramref name="a"/> and <paramref name="b"/>
///		into <paramref name="buffer"/>. </summary>
/// <param name="a"> Start of the edge in normalized device coordinates. </param>
/// <param name="b"> End of the edge in normalized device coordinates. </param>
void Rasterize(const Vec3f& a, const Vec3f& b, BackBuffer& buffer) {
	// This is a really crude algorithm that interpolates between the two points
	// of the edge at a fixed number of points and colors the pixel corresponding
	// to those points.
	constexpr auto steps = 60;
	const auto screenSize = Vec2f(buffer.ColumnCount(), buffer.RowCount());
	for (size_t i = 0; i <= steps; ++i) {
		const float t = float(i) / float(steps);
		const auto ndc = t * a.xy + (1.0f - t) * b.xy; // Interpolated point in NDC.
		const auto screen = (0.5f * ndc + 0.5f) * screenSize; // Interpolated point in screen space.
		const auto pixel = Vec2i(screen + Vec2f(0.5f, 0.5f)); // Round to screen space to integers.

		// Fill pixel if it falls within the bounds of the back buffer.
		if (0 <= pixel.x && pixel.x < buffer.ColumnCount()
			&& 0 <= pixel.y && pixel.y < buffer.RowCount()) {
			buffer(pixel.y, pixel.x) = '#';
		}
	}
}


/// <summary> Rasterize a triangle into <paramref name="buffer"/>. </summary>
void Rasterize(const Triangle3f& triangle, BackBuffer& buffer) {
	// Only the edges are rasterized, the triangle is not filled.
	// A.k.a wireframe rendering.
	if (IsFrontFacing(triangle)) {
		Rasterize(triangle.corners[0], triangle.corners[1], buffer);
		Rasterize(triangle.corners[0], triangle.corners[2], buffer);
		Rasterize(triangle.corners[1], triangle.corners[2], buffer);
	}
}


/// <summary> Render an object into <paramref name="buffer"/>. </summary>
/// <param name="view"> The view transform. </param>
/// <param name="proj"> The projection transform. </param>
void Render(const Object& object, const Mat44f& view, const Mat44f& proj, BackBuffer& buffer) {
	const auto model = object.GetTransform();
	const auto mvp = model * view * proj;
	for (const auto& tri : object.triangles) {
		const auto transformed = TransformTriangle(tri, mvp);
		Rasterize(transformed, buffer);
	}
}


/// <summary> Creates a unit cube. </summary>
Object Cube() {
	// A single face of the unit cube.
	const std::array face = {
		Triangle3f(Vec3f(-1, 1, 1), Vec3f(1, -1, 1), Vec3f(-1, -1, 1)),
		Triangle3f(Vec3f(-1, 1, 1), Vec3f(1, 1, 1), Vec3f(1, -1, 1)),
	};

	// The single face is rotated in 6 different ways to create the 6 faces of the cube.
	std::vector triangles = {
		face[0],
		face[1],
		TransformTriangle(face[0], RotationY(Deg2Rad(90.f))),
		TransformTriangle(face[1], RotationY(Deg2Rad(90.f))),
		TransformTriangle(face[0], RotationY(Deg2Rad(180.f))),
		TransformTriangle(face[1], RotationY(Deg2Rad(180.f))),
		TransformTriangle(face[0], RotationY(Deg2Rad(270.f))),
		TransformTriangle(face[1], RotationY(Deg2Rad(270.f))),
		TransformTriangle(face[0], RotationX(Deg2Rad(90.f))),
		TransformTriangle(face[1], RotationX(Deg2Rad(90.f))),
		TransformTriangle(face[0], RotationX(Deg2Rad(-90.f))),
		TransformTriangle(face[1], RotationX(Deg2Rad(-90.f))),
	};

	// The object's transform is set to identity.
	return Object{ std::move(triangles), { 1, 1, 1 }, Identity(), { 0, 0, 0 } };
}


/// <summary> Clears the back buffer. </summary>
void Clear(BackBuffer& buffer) {
	for (size_t row = 0; row < buffer.RowCount(); ++row) {
		for (size_t col = 0; col < buffer.ColumnCount(); ++col) {
			buffer(row, col) = ' ';
		}
	}
}


/// <summary> Display the back buffer on STDOUT. </summary>
void Display(const BackBuffer& buffer) {
	// Clears the terminal, hopefully.
#ifdef _WIN32
	system("cls");
#elif defined(__linux__)
	system("clear");
#endif

	// Draw a little rim around the rendered buffer.
	std::cout << '+';
	for (size_t col = 0; col < buffer.ColumnCount(); ++col) {
		std::cout << '-';
	}
	std::cout << '+' << '\n';

	for (ptrdiff_t row = buffer.RowCount() - 1; row >= 0; --row) {
		std::cout << '|'; // The left side of the rim.

		// Display the "pixels" of the back buffer.
		for (size_t col = 0; col < buffer.ColumnCount(); ++col) {
			std::cout << buffer(row, col);
		}

		std::cout << '|' << '\n'; // The right side of the rim.
	}

	// The bottom part of the rim.
	std::cout << '+';
	for (size_t col = 0; col < buffer.ColumnCount(); ++col) {
		std::cout << '-';
	}
	std::cout << '+' << '\n';
}


int main() {
	BackBuffer buffer;
	Camera camera;

	// Set up the camera.
	const auto screenAspect = float(buffer.ColumnCount()) / float(buffer.RowCount());
	const auto characterAspect = 0.5f; // Characters are not square, like pixels.
	const auto effectiveAspect = screenAspect * characterAspect;

	camera.SetEye({ 3, 2, 2 });
	camera.SetTarget({ 0, 0, 0 });
	camera.SetUp({ 0, 0, 1 });
	camera.SetAspect(effectiveAspect);
	camera.SetFoVH(57.f);

	// Declare our game objects.
	auto cube = Cube();

	// Render loop runs for a fixed number of frames. "Good enough for a tech demo."
	for (size_t frame = 0; frame < 10000; ++frame) {
		// Update game objects.
		cube.rotation = Quatf(RotationZ(0.02f)) * cube.rotation; // Slowly spin the cube around.

		// Render game objects.
		Clear(buffer);
		const auto view = camera.GetViewTransform();
		const auto proj = camera.GetPerspectiveTransform();
		Render(cube, view, proj, buffer);
		Display(buffer);

		// Maintain about 30 FPS.
		std::this_thread::sleep_for(std::chrono::milliseconds(33));
	}
}