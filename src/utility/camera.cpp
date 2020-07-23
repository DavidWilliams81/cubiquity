#include "camera.h"

using namespace Cubiquity;

Camera::Camera()
	:position(0.0f, 0.0f, 0.0f)
	,pitch(0.0f)
	,yaw(0.0f)
{
	fovInDegrees = 60.0f; // FIXME - Shouldn't be hard-coded.
	aspect = 4.0f / 3.0f; // FIXME - Shouldn't be hard-coded.
}

Ray3d Camera::rayFromViewportPos(int x, int y, int width, int height) const
{
	Vector3d camPos = static_cast<Vector3d>(position);
	Vector3d camDir = static_cast<Vector3d>(forward());
	Vector3d camUp = static_cast<Vector3d>(up());
	Vector3d camRight = static_cast<Vector3d>(right());

	double invWidth = 1.0f / width;
	double invHeight = 1.0f / height;

	float aspectRatio = static_cast<float>(width) / static_cast<float>(height); // FIXME - we already have this stored in the camera.

	float scale = tan(fovInDegrees * 0.0174533f * 0.5f) * 2.0f;

	float xOffset = x - (width / 2.0f) + 0.5f;
	float yOffset = y - (height / 2.0f) + 0.5f;

	Vector3d camTarget = camPos + camDir;
	camTarget += camRight * (invWidth * xOffset * aspectRatio * scale);
	camTarget -= camUp * (invHeight * yOffset * scale); // Inverts Y

	Vector3d rayDir = camTarget - camPos;
	rayDir = normalize(rayDir);
	Ray3d ray(camPos, rayDir);

	return ray;
}

Vector3f Camera::forward() const
{
	// Direction : Spherical coordinates to Cartesian coordinates conversion
	return Vector3f
	(
		cos(pitch) * sin(yaw),
		sin(pitch),
		cos(pitch) * cos(yaw)
	);
}

Vector3f Camera::right() const
{
	// Right vector
	return Vector3f(
		sin(yaw - 3.14f / 2.0f),
		0,
		cos(yaw - 3.14f / 2.0f)
	);
}

Vector3f Camera::up() const
{
	// Up vector
	return cross(right(), forward());
}

Matrix4x4f Camera::viewMatrix() const
{
	// Camera matrix
	return lookAtRH(
		position,           // Camera is here
		position + forward(), // and looks here : at the same position, plus "direction"
		up()                  // Head is up (set to 0,-1,0 to look upside-down)
	);
}

Matrix4x4f Camera::projectionMatrix() const
{
	return perspective_matrix(fovInDegrees * 0.0174533f, aspect, 0.1f, 10000.0f);
}
