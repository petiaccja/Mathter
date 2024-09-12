#pragma once

#include "Math.hpp"


/// <summary> A camera object that help with view and projective transforms. </summary>
class Camera {
public:
	void SetEye(const Vec3f& eye);
	const Vec3f& GetEye() const;

	void SetTarget(const Vec3f& target);
	const Vec3f& GetTarget() const;

	void SetUp(const Vec3f& eye);
	const Vec3f& GetUp() const;

	void SetFoVH(float fovH);
	float GetFoVH() const;

	void SetAspect(float aspect);
	float GetAspect() const;

	Mat44f GetViewTransform() const;
	Mat44f GetPerspectiveTransform() const;

private:
	Vec3f m_eye = { 0, 0, 0 };
	Vec3f m_target = { 1, 0, 0 };
	Vec3f m_up = { 0, 0, 1 };
	float m_fovH = 75.0f;
	float m_aspect = 16.0f / 9.0f;
};