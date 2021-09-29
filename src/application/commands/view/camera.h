#ifndef CAMERA_H_1B80A34E
#define CAMERA_H_1B80A34E

#include "geometry.h"

class Camera
{
public:
	Camera();

	Cubiquity::Vector3f position;
	float pitch;
	float yaw;
	float fovInDegrees;
	float aspect;

	Cubiquity::Ray3d rayFromViewportPos(int x, int y, int width, int height) const;

	Cubiquity::Vector3f forward() const;
	Cubiquity::Vector3f right() const;
	Cubiquity::Vector3f up() const;

	Cubiquity::Matrix4x4f viewMatrix() const;
	Cubiquity::Matrix4x4f projectionMatrix() const;
};

#endif //CAMERA_H_1B80A34E