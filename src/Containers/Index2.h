#ifndef INDEX_2_H
#define INDEX_2_H

#include "Std/CStdInt.h"

namespace asa
{
class Index2
{
public:
    uint x, y;

    Index2()
        : x(0)
        , y(0)
    {
    }

    Index2(uint x, uint y)
    {
        this->x = x;
        this->y = y;
    }

    uint &operator[](uint i) { return i == 0 ? x : y; }
    const uint &operator[](uint i) const { return i == 0 ? x : y; }
};

bool operator==(const Index2 &lhs, const Index2 &rhs);
bool operator!=(const Index2 &lhs, const Index2 &rhs);
};  // namespace asa

#endif
