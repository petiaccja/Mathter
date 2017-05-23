//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)
#include <gtest/gtest.h>

#include "Mathter/Quaternion.hpp"

using namespace mathter;


TEST(Quaternion, Ctor) {
	Quaternion<float> q(1, 2, 3, 4);
	ASSERT_EQ(q.w, 1);
	ASSERT_EQ(q.x, 2);
	ASSERT_EQ(q.y, 3);
	ASSERT_EQ(q.z, 4);
}


TEST(Quaternion, AxisAngle) {
	
}