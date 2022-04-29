#include "Index2.h"

namespace asa
{
bool operator==(const Index2 &lhs, const Index2 &rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y;
}

bool operator!=(const Index2 &lhs, const Index2 &rhs)
{
    return !(lhs == rhs);
}
}  // namespace asa
