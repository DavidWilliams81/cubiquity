#include "camera.h"

using namespace Cubiquity;

Camera::Camera()
	:position({ 0.0, 0.0, 0.0 })
	,pitch(0.0)
	,yaw(0.0)
{
	fovInDegrees = 60.0; // FIXME - Shouldn't be hard-coded.
	aspect = 4.0 / 3.0; // FIXME - Shouldn't be hard-coded.
}

Ray3d Camera::rayFromViewportPos(int x, int y, int width, int height) const
{
	vec3d camPos = static_cast<vec3d>(position);
	vec3d camDir = static_cast<vec3d>(forward());
	vec3d camUp = static_cast<vec3d>(up());
	vec3d camRight = static_cast<vec3d>(right());

	double invWidth = 1.0f / width;
	double invHeight = 1.0f / height;

	float aspectRatio = static_cast<float>(width) / static_cast<float>(height); // FIXME - we already have this stored in the camera.

	float scale = tan(fovInDegrees * 0.0174533f * 0.5f) * 2.0f;

	float xOffset = x - (width / 2.0f) + 0.5f;
	float yOffset = y - (height / 2.0f) + 0.5f;

	vec3d camTarget = camPos + camDir;
	camTarget += camRight * (invWidth * xOffset * aspectRatio * scale);
	camTarget -= camUp * (invHeight * yOffset * scale); // Inverts Y

	vec3d rayDir = camTarget - camPos;
	rayDir = normalize(rayDir);
	Ray3d ray(camPos, rayDir);

	return ray;
}

vec3d Camera::forward() const
{
	// Direction : Spherical coordinates to Cartesian coordinates conversion
	return vec3d
	({
		// Looking along (0,1,0) when pitch and yaw are both zero.
		cos(pitch) * sin(yaw),
		cos(pitch) * cos(yaw),
		sin(pitch)
	});
}

vec3d Camera::right() const
{
	// Looking along (1,0,0) when pitch and yaw are both zero.
	return vec3d({
		sin(yaw + (Pi / 2)),
		cos(yaw + (Pi / 2)),
		0
	});
}

vec3d Camera::up() const
{
	// Up vector
	return cross(right(), forward());
}

Matrix4x4d Camera::viewMatrix() const
{
	// Camera matrix
	return lookAtRH(
		position,           // Camera is here
		position + forward(), // and looks here : at the same position, plus "direction"
		up()                  // Head is up (set to 0,-1,0 to look upside-down)
	);
}

Matrix4x4d Camera::projectionMatrix() const
{
	return perspective_matrix(fovInDegrees * 0.0174533, aspect, 1.0, 10000.0);
}
