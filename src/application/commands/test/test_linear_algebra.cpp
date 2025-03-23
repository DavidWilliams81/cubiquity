#include "test_linear_algebra.h"

#include "geometry.h"

using namespace Cubiquity;
using namespace std;

bool testLinearAlgebra()
{
	vec2i a({ 1,2 });
	vec2i b({ 3,4 });

	a += 1;
	a = a + 1;

	a += b;
	vec2i c = a + b;

	vec2f f = static_cast<vec2f>(c);

	/*vec3f vf(1.5f, 2.5f, 3.5f);
	vec3i vi = static_cast<vec3i>(vf);
	log_info(vi)*/

	/*Matrix4x4f translation = translationMatrix(vec3f(0.0, 0.0, -10));
	Matrix4x4f rotation = rotationMatrix(vec4f(0.3f, 1.0f, 0.0f, 0.0f));
	Matrix4x4f scale = scalingMatrix(vec3f(1.0f, 2.0f, 3.0f));

	auto result = mul(mul(translation, rotation), scale);

	log_info(result)*/

	return true;
}
