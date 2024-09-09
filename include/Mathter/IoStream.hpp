// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Common/TypeTraits.hpp"
#include "Matrix/Matrix.hpp"
#include "Quaternion/Math.hpp"
#include "Quaternion/Quaternion.hpp"
#include "Utility.hpp"
#include "Vector/Vector.hpp"

#include <iostream>
#include <type_traits>

namespace mathter {


namespace impl {

	template <class Char, class CharTraits>
	char ReadWhileSpace(std::basic_istream<Char, CharTraits>& is) {
		char c = Char('\0');
		do {
			if (!is.good()) {
				throw std::runtime_error("input stream failed");
			}
			c = is.peek();
			if (std::isspace(c)) {
				is.get();
			}
		} while (std::isspace(c));
		return c;
	}


	template <class Char, class CharTraits>
	std::basic_string<Char, CharTraits> ReadArray(std::basic_istream<Char, CharTraits>& is) {
		std::basic_string<Char, CharTraits> str;
		size_t depth = 0;
		do {
			const auto c = static_cast<char>(is.get());
			if (!is.good()) {
				throw std::invalid_argument("unexpected end of stream");
			}
			if (c == Char('[')) {
				depth += 1;
			}
			else if (c == Char(']')) {
				depth -= 1;
			}
			str += c;
		} while (depth > 0);
		return str;
	}


	template <class Char, class CharTraits>
	std::vector<std::string_view> SplitArray(std::basic_string_view<Char, CharTraits> str, Char sep = ',') {
		std::vector<std::basic_string_view<Char, CharTraits>> substrings;
		auto prev = 0;
		size_t depth = 0;
		for (size_t i = 0; i < str.size(); ++i) {
			if (str[i] == '[') {
				depth += 1;
			}
			else if (str[i] == ']') {
				depth -= 1;
			}
			else if (depth == 0 && str[i] == sep) {
				substrings.push_back(str.substr(prev, i - prev));
				prev = i + 1;
			}
		}
		substrings.push_back(str.substr(prev, str.npos));
		return substrings;
	}


	template <class Char, class CharTraits,
			  class Element>
	void PrintElement(std::basic_ostream<Char, CharTraits>& os, const Element& value) {
		if constexpr (std::is_integral_v<Element>) {
			os << value;
		}
		else if constexpr (std::is_floating_point_v<Element>) {
			os << value;
		}
		else if constexpr (is_complex_v<Element>) {
			os << std::real(value);
			if (std::imag(value) != remove_complex_t<Element>(0)) {
				os << " + " << std::imag(value) << "j";
			}
		}
		else {
			os << value;
		}
	}


	template <class Char, class CharTraits,
			  class Element>
	void ParseElement(std::basic_istream<Char, CharTraits>& is, Element& value) {
		if constexpr (std::is_integral_v<Element>) {
			is >> value;
		}
		else if constexpr (std::is_floating_point_v<Element>) {
			is >> value;
		}
		else if constexpr (is_complex_v<Element>) {
			remove_complex_t<Element> re(0);
			remove_complex_t<Element> im(0);
			is >> re;
			const char next = ReadWhileSpace(is);
			if (next == '+') {
				is.get();
				is >> im;
			}
			value = { re, im };
		}
		else {
			is >> value;
		}
	}


	template <class Char, class CharTraits,
			  class Element, size_t Size>
	void PrintArray(std::basic_ostream<Char, CharTraits>& os, const std::array<Element, Size>& array) {
		os << "[";
		for (const auto& element : array) {
			PrintElement(os, element);
			&element != &array.back() ? static_cast<void>(os << ", ") : void();
		}
		os << "]";
	}


	template <class Char, class CharTraits,
			  class Element, size_t Size>
	void ParseArray(std::basic_istream<Char, CharTraits>& is, std::array<Element, Size>& array) {
		const char next = ReadWhileSpace(is);
		if (next != '[') {
			throw std::invalid_argument("expectd '['");
		}

		const auto str = ReadArray(is);

		const auto elements = SplitArray(std::basic_string_view<Char, CharTraits>(str).substr(1, str.size() - 2));
		if (elements.size() != array.size()) {
			throw std::invalid_argument("incorrect number of array elements");
		}
		for (size_t i = 0; i < array.size(); ++i) {
			std::basic_istringstream<Char, CharTraits> iss;
			iss.str(std::basic_string<Char, CharTraits>(elements[i]));
			ParseElement(iss, array[i]);
		}
	}

} // namespace impl


template <class Char, class CharTraits,
		  class T, int Dim, bool Packed>
auto operator<<(std::basic_ostream<Char, CharTraits>& os, const Vector<T, Dim, Packed>& v) -> std::basic_ostream<Char, CharTraits>& {
	const auto packed = Vector<T, Dim, true>(v);
	impl::PrintArray(os, packed.elements.array);
	return os;
}


template <class Char, class CharTraits,
		  class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto operator<<(std::basic_ostream<Char, CharTraits>& os, const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) -> std::basic_ostream<Char, CharTraits>& {
	std::array<Vector<T, Columns, Packed>, Rows> rows;
	for (size_t row = 0; row < Rows; ++row) {
		rows[row] = m.Row(row);
	}
	impl::PrintArray(os, rows);
	return os;
}


template <class Char, class CharTraits,
		  class T, eQuaternionLayout Layout, bool Packed>
auto operator<<(std::basic_ostream<Char, CharTraits>& os, const Quaternion<T, Layout, Packed>& q) -> std::basic_ostream<Char, CharTraits>& {
	os << Char('[')
	   << Rad2Deg(Angle(q)) << " deg"
	   << Char(' ') << Char('@') << Char(' ')
	   << Axis(q)
	   << Char(']');
	return os;
}


template <class Char, class CharTraits,
		  class T, int Dim, bool Packed>
auto operator>>(std::basic_istream<Char, CharTraits>& is, Vector<T, Dim, Packed>& v) -> std::basic_istream<Char, CharTraits>& {
	auto packed = Vector<T, Dim, true>();
	impl::ParseArray(is, packed.elements.array);
	v = Vector<T, Dim, Packed>(packed);
	return is;
}


template <class Char, class CharTraits,
		  class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto operator>>(std::basic_istream<Char, CharTraits>& is, Matrix<T, Rows, Columns, Order, Layout, Packed>& m) -> std::basic_istream<Char, CharTraits>& {
	std::array<Vector<T, Columns, Packed>, Rows> rows;
	impl::ParseArray(is, rows);
	for (size_t row = 0; row < Rows; ++row) {
		m.Row(row, rows[row]);
	}
	return is;
}


template <class Char, class CharTraits,
		  class T, eQuaternionLayout Layout, bool Packed>
auto operator>>(std::basic_istream<Char, CharTraits>& is, Quaternion<T, Layout, Packed>& q) -> std::basic_istream<Char, CharTraits>& {
	const auto str = impl::ReadArray(is);
	const auto elements = impl::SplitArray(std::basic_string_view<Char, CharTraits>(str).substr(1, str.size() - 2), Char('@'));
	if (elements.size() != 2) {
		throw std::invalid_argument("expected an angle and an axis");
	}
	T angle;
	Vector<T, 3, Packed> axis;

	std::basic_stringstream<Char, CharTraits> ss;
	ss.str(std::basic_string<Char, CharTraits>(elements[0]));
	impl::ParseElement(ss, angle);
	ss.str(std::basic_string<Char, CharTraits>(elements[1]));
	ss >> axis;
	const auto theta = Deg2Rad(angle) / T(2);
	q = Quaternion<T, Layout, Packed>(std::cos(theta), axis * std::sin(theta));
	return is;
}

} // namespace mathter