#include "Object.hpp"

#include "Mathter/Transforms.hpp"


Mat44f Object::GetTransform() const {
	// The scaling and translation transforms must be converted to Mat44 before they can be used for multiplication.
	// Since we defined Mat44 to have FOLLOW_VECTOR multiplication order, we're concatenating
	// the transforms from left to right.
	return Mat44f(Scale(scale)) * Mat44f(rotation) * Mat44f(Translation(position));
}