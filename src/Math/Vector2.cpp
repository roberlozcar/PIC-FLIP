#include "Vector2.h"

namespace asa
{
const Vector2 Vector2::ZERO(0, 0);

Vector2 operator*(const float &k, const Vector2 &v)
{
    return Vector2(k * v.x, k * v.y);
}

const Matrix2 Matrix2::ZERO(0, 0, 0, 0);

const Matrix2 Matrix2::IDENTITY(1, 0, 0, 1);

Matrix2 operator*(const float &k, const Matrix2 &m)
{
    return Matrix2(m.v[0][0] * k, m.v[0][1] * k, m.v[1][0] * k, m.v[1][1] * k);
}
}  // namespace asa
