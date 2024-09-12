#pragma once

#include "Math.hpp"

#include <vector>


/// <summary> A polygonal 3D geometry consisting of triangles. </summary>
struct Object {
	std::vector<Triangle3f> triangles;

	Vec3f scale;
	Quatf rotation;
	Vec3f position;

	Mat44f GetTransform() const;
};