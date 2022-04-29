#ifndef VECTOR_3_H
#define VECTOR_3_H

#include "Std/CMath.h"

namespace asa
{
class Vector3
{
public:
    float x, y, z;

public:
    Vector3()
        : x(0)
        , y(0)
        , z(0)
    {
    }

    Vector3(const float x, const float y, const float z)
        : x(x)
        , y(y)
        , z(z)
    {
    }

    void operator*=(const float k)
    {
        x *= k;
        y *= k;
        z *= k;
    }

    void operator+=(const Vector3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    void operator-=(const Vector3 &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    void operator*=(const Vector3 &v)
    {
        x *= v.x;
        y *= v.y;
        z *= v.z;
    }

    Vector3 operator+(const float &k) const { return Vector3(x + k, y + k, z + k); }

    Vector3 operator-(const float &k) const { return Vector3(x - k, y - k, z - k); }

    Vector3 operator*(const float &k) const { return Vector3(x * k, y * k, z * k); }

    Vector3 operator+(const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }

    Vector3 operator-(const Vector3 &v) const { return Vector3(x - v.x, y - v.y, z - v.z); }

    Vector3 operator*(const Vector3 &v) const { return Vector3(x * v.x, y * v.y, z * v.z); }

    Vector3 operator/(const Vector3 &v) const { return Vector3(x / v.x, y / v.y, z / v.z); }

    float dot(const Vector3 &v) const { return x * v.x + y * v.y + z * v.z; }

    float lengthSqr() const { return x * x + y * y + z * z; }

    float length() const { return std::sqrt(lengthSqr()); }

public:
    const static Vector3 ZERO;
};

template <>
inline Vector3 zeroConstant<Vector3>()
{
    return Vector3::ZERO;
}
}  // namespace asa

#endif
