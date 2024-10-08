# Changelog

## 2.0.0

- Version 2.0.0 is a complete overhaul:
	- The interfaces stayed mostly the same (especially the common ones)
	- The internals have been completely rewriten and/or refactored
	- Migrating from 1.x.x:
		- Read this changelog
		- Spam recompile until it works
- Folder structure and library organization:
	- Changes to the top-level headers:
		- `<Mathter/Matrix.hpp>` no longer includes decompositions (include decompositions individually)
		- `<Mathter/Matrix.hpp>` no longer includes transforms (include `<Mathter/Transforms.hpp>`)
	- The internals (e.g. `<Mathter/Vector/*.hpp>`) have been completely refactored
- Vectors:
	- New data structure:
		- No more union abuse, `reinterpret_cast`, and UB
		- Named accessors (e.g. `.x`) are no longer T&, but are implemented as swizzlers
	- Homogeneous downcast now performs perspective division
	- CTAD support added to some constructors
	- Deprecated concatenation by the pipe (`|`) operator
		- Use concatenating constructor with CTAD instead
	- `IsNullvector`: removed due to questionable meaning / use-cases
	- `IsNormalized`: removed due to questionable meaning / use-cases
	- `SafeNormalize`: renamed to `NormalizePrecise`
- Matrices:
	- Added additional converting constructors (switch order, layout, and scalar type)
	- Removed casting utilities:
		- Use the converting constructors instead
		- Similar utilities in `Cast.hpp` (warning: not intended for end users)
	- Removed submatrices / matrix views
	- Added `Row` and `Column` accessors (these are layout-agnostic)
	- Added `Extract` and `Insert` accessors (these replace sub-matrices)
- Quaternions:
	- Added template parameter to select memory layout:
		- Choose between SIJK and IJKS ordering
	- New data structure:
		- Deprecated `VectorPart` and `ScalarPart`
			- Use the swizzlers `.scalar` and `.vector`
		- Named accessors (e.g. `.s`) are no longer T&, but are implemented as swizzlers
	- `Axis` and `Angle` are now free functions
	- Renamed `Conjugate` to `Conj` (it's different from the quaternion inverse)
	- `Inverse` does actual inverse (same as conjugation **only** for unit quaternions)
	- Removed conversion to `Vector<T, 3>` (ambiguous semantics)
	- Deprecated `operator*(Quaternion, Vector)`:
		- Use `Quaternion::operator()(Vector)` to apply rotations to vectors
- Transforms:
	- New transforms:
		- Generalized shear
		- Random (this is not really a transform though)
	- Transforms applicable to more types (e.g. `Vector<float, 3> v = Zero()` is now valid)
- Geometry:
	- The intersections have been completely reworked
		- Use `Intersect(a, b)`: returns an `optional` of the intersection's coordinate
- Decompositions:
	- Unified interface for the different decompositions:
		- `Solve`: solve one/many equation systems & least squares problems
		- `Inverse`: compute inverse and pseudoinverse
		- Get factors via structured bindings
	- QR:
		- Added LQ decomposition
	- SVD:
		- Added 1-sided Jacobi algorithm (alongside existing 2-sided)
- IoStream:
	- Complete rewrite (old implementation was mostly a quick and dirty utility)
	- Changed formatting to Python's array syntax
- General:
	- Mixed precision now widely permitted (e.g. adding vector or doubles to vector of complex floats returns a vector of complex doubles)
	- Full support for complex numbers:
		- Vectors and matrices now fully support complex numbers
		- Methods (like `Dot`) have been adapted
		- Algorithms (like decompositions) have been adapted
	- Improved testing:
		- More accurately (and conservatively) measuring coverage
		- Raised coverage above 90%, with conservative estimates
		- More template parameter combinations are tested (e.g. complex numbers, mixed matrix layouts)
	- Many new mathematical functions: most on vectors and matrices, a few on quaternions
- Performance:
	- Added many new benchmarks
- Documentation:
	- Added an examples folder with a mini project
	- Rewrote the guide