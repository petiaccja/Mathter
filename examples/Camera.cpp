#include "Camera.hpp"

#include <Mathter/Transforms.hpp>
#include <Mathter/Utility.hpp>


void Camera::SetEye(const Vec3f& eye) {
	m_eye = eye;
}


const Vec3f& Camera::GetEye() const {
	return m_eye;
}


void Camera::SetTarget(const Vec3f& target) {
	m_target = target;
}


const Vec3f& Camera::GetTarget() const {
	return m_target;
}


void Camera::SetUp(const Vec3f& up) {
	m_up = up;
}


const Vec3f& Camera::GetUp() const {
	return m_up;
}


void Camera::SetFoVH(float fovH) {
	m_fovH = fovH;
}


float Camera::GetFoVH() const {
	return m_fovH;
}


void Camera::SetAspect(float aspect) {
	m_aspect = aspect;
}


float Camera::GetAspect() const {
	return m_aspect;
}


Mat44f Camera::GetViewTransform() const {
	// This creates a left-handed camera space where the Z axis is forward,
	// the X axis if right, and the Y axis is up.
	return LookAt(m_eye, m_target, m_up, true, false, false);
}


Mat44f Camera::GetPerspectiveTransform() const {
	// This projects the camera space to normalized device coordinates.
	// If the camera space used negative Z forward, the near and far planes should also be negative.
	// The normalized device coordinates are X in [-1, 1], Y in [-1, 1], and Z in [0, 1].
	return mathter::Perspective(mathter::Deg2Rad(m_fovH), m_aspect, 0.1f, 100.f, 0.0f, 1.0f);
}
