#ifndef CAMERA_H_1B80A34E
#define CAMERA_H_1B80A34E

#include "base/ray.h"

const float Pi = 3.14159265358979f;

class Camera
{
public:
	Camera();

	dvec3 position;
	double pitch;
	double yaw;
	double fovInDegrees;
	double aspect;

	Ray3d rayFromViewportPos(int x, int y, int width, int height) const;

	dvec3 forward() const;
	dvec3 right() const;
	dvec3 up() const;

	dmat4 viewMatrix() const;
	dmat4 projectionMatrix() const;
};

#endif //CAMERA_H_1B80A34E