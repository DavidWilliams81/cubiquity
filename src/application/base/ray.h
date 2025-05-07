#ifndef CUBIQUITY_APP_RAY_H
#define CUBIQUITY_APP_RAY_H

#include "base/types.h"

template <class VecType>
class Ray3
{
public:
	Ray3() {}

	Ray3(const VecType& origin, const VecType& dir)
		: mOrigin(origin)
		, mDir(dir) {
	}

	template <typename CastType> explicit Ray3(const Ray3<CastType>& ray)
	{
		//for (u32 ct = 0; ct < Size; ++ct) { data[ct] = static_cast<Type>(vector.data[ct]); }
		mOrigin = static_cast<VecType>(ray.mOrigin);
		mDir = static_cast<VecType>(ray.mDir);
	}

public:
	VecType mOrigin;
	VecType mDir; // FIXME - Would float always be sufficient for the direction?
};

typedef Ray3<vec3> Ray3f;
typedef Ray3<dvec3> Ray3d;

#endif // CUBIQUITY_APP_RAY_H
