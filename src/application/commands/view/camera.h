#ifndef CAMERA_H_1B80A34E
#define CAMERA_H_1B80A34E

#include "geometry.h"

class Camera
{
public:
	Camera();

	Cubiquity::Vector3d position;
	double pitch;
	double yaw;
	double fovInDegrees;
	double aspect;

	Cubiquity::Ray3d rayFromViewportPos(int x, int y, int width, int height) const;

	Cubiquity::Vector3d forward() const;
	Cubiquity::Vector3d right() const;
	Cubiquity::Vector3d up() const;

	Cubiquity::Matrix4x4d viewMatrix() const;
	Cubiquity::Matrix4x4d projectionMatrix() const;
};

#endif //CAMERA_H_1B80A34E