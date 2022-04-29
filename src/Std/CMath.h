#ifndef ASA_CMATH_H
#define ASA_CMATH_H

#include <Std/CStdInt.h>
#include <cassert>
#include <cmath>

namespace asa
{
template <typename T>
T lerp(const T &a, const T &b, const float t)
{
    assert(0 <= t && t <= 1);

    return a * (1.0f - t) + b * t;
}

template <typename T>
T bilerp(const T &aa, const T &ba, const T &ab, const T &bb, const float tx, const float ty)
{
    const T y1 = lerp(aa, ba, tx);
    const T y2 = lerp(ab, bb, tx);
    return lerp(y1, y2, ty);
}

template <typename T>
T clamp(const T &v, const T &lo, const T &hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

inline uint clamp(const uint &v, const uint &lo, const uint &hi)
{
    const int res = clamp(static_cast<int>(v), static_cast<int>(lo), static_cast<int>(hi));
    return static_cast<uint>(res);
}

template <typename T>
T zeroConstant()
{
    return static_cast<T>(0);
}
}  // namespace asa

#endif
