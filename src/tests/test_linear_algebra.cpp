#include "test_linear_algebra.h"

#include "geometry.h"

using namespace Cubiquity;
using namespace std;

bool testLinearAlgebra()
{
	Vector2i a(1,2);
	Vector2i b(3,4);

	a += 1;
	a = a + 1;

	a += b;
	Vector2i c = a + b;

	/*Vector3f vf(1.5f, 2.5f, 3.5f);
	Vector3i vi = static_cast<Vector3i>(vf);
	std::cout << vi << std::endl;*/

	/*Matrix4x4f translation = translationMatrix(Vector3f(0.0, 0.0, -10));
	Matrix4x4f rotation = rotationMatrix(Vector4f(0.3f, 1.0f, 0.0f, 0.0f));
	Matrix4x4f scale = scalingMatrix(Vector3f(1.0f, 2.0f, 3.0f));

	auto result = mul(mul(translation, rotation), scale);

	std::cout << std::fixed;
	std::cout << std::setprecision(3);
	std::cout << result << std::endl;*/

	return true;
}
