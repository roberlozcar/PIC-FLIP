#ifndef VECTOR_2_H
#define VECTOR_2_H

#include "Std/CMath.h"

namespace asa
{
class Vector2
{
public:
    float x, y;

public:
    Vector2()
        : x(0.0f)
        , y(0.0f)
    {
    }

    Vector2(const float x_, const float y_)
        : x(x_)
        , y(y_)
    {
    }

    void operator+=(const Vector2 &v)
    {
        x += v.x;
        y += v.y;
    }

    void operator-=(const Vector2 &v)
    {
        x -= v.x;
        y -= v.y;
    }

    void operator*=(const Vector2 &v)
    {
        x *= v.x;
        y *= v.y;
    }

    Vector2 operator+(const float &k) const { return Vector2(x + k, y + k); }

    Vector2 operator-(const float &k) const { return Vector2(x - k, y - k); }

    Vector2 operator*(const float &k) const { return Vector2(x * k, y * k); }

    Vector2 operator+(const Vector2 &v) const { return Vector2(x + v.x, y + v.y); }

    Vector2 operator-(const Vector2 &v) const { return Vector2(x - v.x, y - v.y); }

    Vector2 operator*(const Vector2 &v) const { return Vector2(x * v.x, y * v.y); }

    Vector2 operator/(const Vector2 &v) const { return Vector2(x / v.x, y / v.y); }

    float dot(const Vector2 &v) const { return x * v.x + y * v.y; }

    float cross(const Vector2 &v) const { return x * v.y - y * v.x; }

    Vector2 rotate90() const { return Vector2(-y, x); }

    float lengthSqr() const { return x * x + y * y; }

    float length() const { return std::sqrt(lengthSqr()); }

public:
    const static Vector2 ZERO;

    static Vector2 crossProd(const Vector2 &v) { return Vector2(-v.y, v.x); }
};

Vector2 operator*(const float &k, const Vector2 &v);

class Matrix2
{
public:
    float v[2][2];

    Matrix2()
    {
        v[0][0] = 0.0f;
        v[0][1] = 0.0f;
        v[1][0] = 0.0f;
        v[1][1] = 0.0f;
    }
    Matrix2(float a, float b, float c, float d)
    {
        v[0][0] = a;
        v[0][1] = b;
        v[1][0] = c;
        v[1][1] = d;
    }

    void operator+=(const Matrix2 &m)
    {
        v[0][0] += m.v[0][0];
        v[0][1] += m.v[0][1];
        v[1][0] += m.v[1][0];
        v[1][1] += m.v[1][1];
    }
    void operator-=(const Matrix2 &m)
    {
        v[0][0] -= m.v[0][0];
        v[0][1] -= m.v[0][1];
        v[1][0] -= m.v[1][0];
        v[1][1] -= m.v[1][1];
    }

    Matrix2 operator+(const Matrix2 &m) const
    {
        return Matrix2(v[0][0] + m.v[0][0], v[0][1] + m.v[0][1], v[1][0] + m.v[1][0], v[1][1] + m.v[1][1]);
    }

    Matrix2 operator-(const Matrix2 &m) const
    {
        return Matrix2(v[0][0] - m.v[0][0], v[0][1] - m.v[0][1], v[1][0] - m.v[1][0], v[1][1] - m.v[1][1]);
    }

    Vector2 operator*(const Vector2 &m) const
    {
        return Vector2(v[0][0] * m.x + v[0][1] * m.y, v[1][0] * m.x + v[1][1] * m.y);
    }

    Matrix2 operator*(const float &k) const { return Matrix2(v[0][0] * k, v[0][1] * k, v[1][0] * k, v[1][1] * k); }

    Matrix2 operator*(const Matrix2 &m) const
    {
        return Matrix2(v[0][0] * m.v[0][0] + v[0][1] * m.v[1][0],
                       v[0][0] * m.v[0][1] + v[0][1] * m.v[1][1],
                       v[1][0] * m.v[0][0] + v[1][1] * m.v[1][0],
                       v[1][0] * m.v[0][1] + v[1][1] * m.v[1][1]);
    }

    const static Matrix2 ZERO;

    const static Matrix2 IDENTITY;

    friend Matrix2 operator*(float k, const Matrix2 &m)
    {
        return Matrix2(k * m.v[0][0], k * m.v[0][1], k * m.v[1][0], k * m.v[1][1]);
    }
};

Matrix2 operator*(const float &k, const Matrix2 &m);

template <>
inline Vector2 zeroConstant<Vector2>()
{
    return Vector2::ZERO;
}
}  // namespace asa

#endif
