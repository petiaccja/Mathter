//==============================================================================
// This software is distributed under The Unlicense. 
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable: 4244)

#include <Catch2/catch.hpp>

#include "Mathter/Geometry.hpp"


TEST_CASE("Bezier curve", "[Geometry]") {
    BezierCurve<float, 2, 2> curve = {
        Vector<float, 2>{},
        Vector<float, 2>{},
        Vector<float, 2>{}
    };


}