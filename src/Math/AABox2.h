#ifndef A_A_BOX_2_H
#define A_A_BOX_2_H

#include "Vector2.h"

namespace asa
{
class AABox2
{
public:
    Vector2 minPosition;
    Vector2 maxPosition;

    AABox2()
        : minPosition(0, 0)
        , maxPosition(0, 0)
    {
    }

    AABox2(const float &minx, const float &miny, const float &maxx, const float &maxy)
    {
        minPosition = Vector2(minx, miny);
        maxPosition = Vector2(maxx, maxy);
    }

    AABox2(const Vector2 &minpos, const Vector2 &maxpos)
    {
        minPosition = minpos;
        maxPosition = maxpos;
    }
};
}  // namespace asa

#endif
