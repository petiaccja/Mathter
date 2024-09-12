#pragma once

#include <Mathter/Geometry/Triangle.hpp>
#include <Mathter/Matrix.hpp>
#include <Mathter/Quaternion.hpp>
#include <Mathter/Vector.hpp>


// Type aliases for Mathter's types, defining the Example project's conventions.
using Vec2i = mathter::Vector<int, 2, false>;
using Vec2f = mathter::Vector<float, 2, false>;
using Vec3f = mathter::Vector<float, 3, false>;
using Mat44f = mathter::Matrix<float, 4, 4, mathter::eMatrixOrder::FOLLOW_VECTOR, mathter::eMatrixLayout::ROW_MAJOR, false>;
using Quatf = mathter::Quaternion<float, mathter::eQuaternionLayout::SCALAR_FIRST, false>;

using Triangle3f = mathter::Triangle<float, 3>;