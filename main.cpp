#include "Mathter\Vector.hpp"

#include <iostream>
using namespace std;



int main() {
	Vector<float, 3> v(1.0f, 2.0f, 3.0f);
	Vector<float, 3> u(1);
	Vector<float, 4> w(v, 1.2f);
	Vector<float, 3> c(v);
	v.Set(1, 2, 3);
	w.Set(2.5f, v);

	auto d = Vector<float, 3>::Cross(v, u);
	
	return 0;
}