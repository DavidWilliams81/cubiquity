#include "camera.h"

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
	dvec3 camPos = static_cast<dvec3>(position);
	dvec3 camDir = static_cast<dvec3>(forward());
	dvec3 camUp = static_cast<dvec3>(up());
	dvec3 camRight = static_cast<dvec3>(right());

	double invWidth = 1.0f / width;
	double invHeight = 1.0f / height;

	float aspectRatio = static_cast<float>(width) / static_cast<float>(height); // FIXME - we already have this stored in the camera.

	float scale = tan(fovInDegrees * 0.0174533f * 0.5f) * 2.0f;

	float xOffset = x - (width / 2.0f) + 0.5f;
	float yOffset = y - (height / 2.0f) + 0.5f;

	dvec3 camTarget = camPos + camDir;
	camTarget += camRight * (invWidth * xOffset * aspectRatio * scale);
	camTarget -= camUp * (invHeight * yOffset * scale); // Inverts Y

	dvec3 rayDir = camTarget - camPos;
	rayDir = normalize(rayDir);
	Ray3d ray(camPos, rayDir);

	return ray;
}

dvec3 Camera::forward() const
{
	// Direction : Spherical coordinates to Cartesian coordinates conversion
	return dvec3
	({
		// Looking along (0,1,0) when pitch and yaw are both zero.
		cos(pitch) * sin(yaw),
		cos(pitch) * cos(yaw),
		sin(pitch)
	});
}

dvec3 Camera::right() const
{
	// Looking along (1,0,0) when pitch and yaw are both zero.
	return dvec3({
		sin(yaw + (Pi / 2)),
		cos(yaw + (Pi / 2)),
		0
	});
}

dvec3 Camera::up() const
{
	// Up vector
	return cross(right(), forward());
}

dmat4 Camera::viewMatrix() const
{
	// Camera matrix
	return lookat_matrix(
		position,           // Camera is here
		position + forward(), // and looks here : at the same position, plus "direction"
		up()                  // Head is up (set to 0,-1,0 to look upside-down)
	);
}

dmat4 Camera::projectionMatrix() const
{
	return linalg::perspective_matrix(fovInDegrees * 0.0174533, aspect, 1.0, 10000.0);
}
