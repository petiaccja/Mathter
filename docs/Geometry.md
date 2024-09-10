# Geometric primitives

Mathter includes the following geometric primitives:
- Lines
- Line segments
- Rays
- Planes (and hyperplanes)
- Triangles
- Bezier curves

The implementations are rather simple because the main use is to calculate intersections. Mathter includes the `Intersect` function that has several overloads for different pairs of geometric primitives.

This module is intended for implementing support features, such as intersecting GUI elements or in-game objects. It's probably a bad idea to base a full-fledged raytracer on these algorithms.
