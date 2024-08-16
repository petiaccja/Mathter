// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Vector/Vector.hpp"

#include <iostream>

namespace mathter {

enum class eEnclosingBracket {
	NONE,
	PARANTHESE,
	BRACKET,
	BRACE,
};

/// <summary> Prints the vector like [1,2,3]. </summary>
template <class T, int Dim, bool Packed>
std::ostream& operator<<(std::ostream& os, const mathter::Vector<T, Dim, Packed>& v) {
	os << "[";
	for (int x = 0; x < Dim; ++x) {
		os << v(x) << (x == Dim - 1 ? "" : ", ");
	}
	os << "]";
	return os;
}


namespace impl {
	template <class T>
	struct dependent_false {
		static constexpr bool value = false;
	};
	template <class T>
	constexpr bool dependent_false_v = dependent_false<T>::value;

	template <class AritT, typename std::enable_if<std::is_integral<AritT>::value && std::is_signed<AritT>::value, int>::type = 0>
	AritT strtonum(const char* str, const char** end) {
		AritT value;
		value = (AritT)strtoll(str, (char**)end, 10);
		return value;
	}
	template <class AritT, typename std::enable_if<std::is_integral<AritT>::value && !std::is_signed<AritT>::value, int>::type = 0>
	AritT strtonum(const char* str, const char** end) {
		AritT value;
		value = (AritT)strtoull(str, (char**)end, 10);
		return value;
	}
	template <class AritT, typename std::enable_if<std::is_floating_point<AritT>::value, int>::type = 0>
	AritT strtonum(const char* str, const char** end) {
		AritT value;
		value = (AritT)strtold(str, (char**)end);
		return value;
	}

	inline const char* StripSpaces(const char* str) {
		while (*str != '\0' && isspace(*str))
			++str;
		return str;
	};

} // namespace impl

/// <summary> Parses a vector from a string. </summary>
template <class Vec>
Vec strtovec(const char* str, const char** end) {
	Vec ret;

	const char* strproc = str;

	// parse initial bracket if any
	strproc = impl::StripSpaces(strproc);
	if (*strproc == '\0') {
		*end = str;
		return ret;
	}

	char startBracket = *strproc;
	char endBracket;
	bool hasBrackets = false;
	switch (startBracket) {
		case '(':
			endBracket = ')';
			hasBrackets = true;
			++strproc;
			break;
		case '[':
			endBracket = ']';
			hasBrackets = true;
			++strproc;
			break;
		case '{':
			endBracket = '}';
			hasBrackets = true;
			++strproc;
			break;
	}

	// parse elements
	for (int i = 0; i < dimension_v<Vec>; ++i) {
		const char* elemend;
		const auto elem = impl::strtonum<scalar_type_t<Vec>>(strproc, &elemend);
		if (elemend == strproc) {
			*end = str;
			return ret;
		}
		else {
			ret[i] = elem;
			strproc = elemend;
		}
		strproc = impl::StripSpaces(strproc);
		if (*strproc == ',') {
			++strproc;
		}
	}

	// parse ending bracket corresponding to initial bracket
	if (hasBrackets) {
		strproc = impl::StripSpaces(strproc);
		if (*strproc != endBracket) {
			*end = str;
			return ret;
		}
		++strproc;
	}

	*end = strproc;
	return ret;
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) {
	os << "[";
	for (int i = 0; i < mat.Height(); ++i) {
		for (int j = 0; j < mat.Width(); ++j) {
			os << mat(i, j) << (j == mat.Width() - 1 ? "" : ", ");
		}
		if (i < Rows - 1) {
			os << "; ";
		}
	}
	os << "]";
	return os;
}


template <class Mat>
Mat strtomat(const char* str, const char** end) {
	using Vec = Vector<scalar_type_t<Mat>, column_count_v<Mat>, is_packed_v<Mat>>;
	Mat ret;

	const char* strproc = str;

	// parse initial bracket if any
	strproc = impl::StripSpaces(strproc);
	if (*strproc == '\0') {
		*end = str;
		return ret;
	}

	char startBracket = *strproc;
	char endBracket;
	bool hasBrackets = false;
	switch (startBracket) {
		case '(':
			endBracket = ')';
			hasBrackets = true;
			++strproc;
			break;
		case '[':
			endBracket = ']';
			hasBrackets = true;
			++strproc;
			break;
		case '{':
			endBracket = '}';
			hasBrackets = true;
			++strproc;
			break;
	}

	// parse rows
	for (int i = 0; i < row_count_v<Mat>; ++i) {
		const char* rowend;
		const auto row = strtovec<Vec>(strproc, &rowend);
		if (rowend == strproc) {
			*end = str;
			return ret;
		}
		else {
			ret.Row(i) = row;
			strproc = rowend;
		}
		strproc = impl::StripSpaces(strproc);
		if (i < row_count_v<Mat> - 1) {
			if (*strproc == ';') {
				++strproc;
			}
			else {
				*end = str;
				return ret;
			}
		}
	}

	// parse ending bracket corresponding to initial bracket
	if (hasBrackets) {
		strproc = impl::StripSpaces(strproc);
		if (*strproc != endBracket) {
			*end = str;
			return ret;
		}
		++strproc;
	}

	*end = strproc;
	return ret;
}


template <class T, eQuaternionLayout Layout, bool Packed>
std::ostream& operator<<(std::ostream& os, const Quaternion<T, Layout, Packed>& q) {
	os << "["
	   << Angle(q) * T(180.0) / T(3.1415926535897932384626)
	   << " deg @ "
	   << Axis(q)
	   << "]";
	return os;
}


} // namespace mathter