### Version 2.0 (?)

Changes:
- Refactored utility folder structure and code (`Common/*`)
- Refactored vector folder structure and code (`Vector/*`)
- `Vector`s:
	- Memory representation:
		- Removed undefined behavior
		- Rewored `Swizzle`s
		- Named accessors (i.e. `.x`) are now `Swizzle`s instead of `T&`
	- CTAD for `Vector` constructors
	- Better support for heterogeneous types (i.e. adding vector of `double` plus vector of `float`)
	- Better support and testing for `std::complex`
	- Functions:
		- `IsNullvector`: removed due to questionable meaning / use-cases
		- `IsNormalized`: removed due to questionable meaning / use-cases
		- `DistancePrecise`: new function
		- `SafeNormalize`: renamed to `NormalizePrecise`
		- Homogeneous downcast now performs perspective division