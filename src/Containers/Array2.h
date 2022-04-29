#ifndef ARRAY_2_H
#define ARRAY_2_H

#include "Index2.h"

namespace asa
{
template <class T>
class Array2
{
public:
    Array2()
        : ptr(nullptr)
    {
    }

    Array2(const uint sx, const uint sy)
        : ptr(nullptr)
    {
        resize(Index2(sx, sy));
    }

        Array2(const Index2 &size)
        : ptr(nullptr)
    {
        resize(size);
    }

    Array2(const Array2<T> &array)
        : ptr(nullptr)
    {
        copy(array);
    }

    ~Array2()
    {
        if (ptr)
            delete[] ptr;
    }

    const Index2 &getSize() const { return size; }

    const T *getData() const { return ptr; }

    uint getLinearIndex(const uint i, const uint j) const { return j * size.x + i; }

    bool resize(const Index2 &size_)
    {
        if (size != size_) {
            if (ptr)
                delete[] ptr;

            size = size_;
            ptr = new T[size.x * size.y];

            clear();

            return true;
        }
        return false;
    }

    void copy(const Array2<T> &src)
    {
        resize(src.size);

        for (uint i = 0, n = size.x * size.y; i < n; ++i)
            ptr[i] = src[i];
    }

    void operator=(const Array2<T> &src) { copy(src); }

    const T &getValue(const Index2 &id) const
    {
        const uint idx = getLinearIndex(id.x, id.y);
        return ptr[idx];
    }
    const T &getValue(const uint i, const uint j) const
    {
        const uint idx = getLinearIndex(i, j);
        return ptr[idx];
    }

    void setValue(const Index2 &id, const T &value)
    {
        const uint idx = getLinearIndex(id.x, id.y);
        ptr[idx] = value;
    }
    void setValue(const uint i, const uint j, const T &value)
    {
        const uint idx = getLinearIndex(i, j);
        ptr[idx] = value;
    }

    const T &operator[](uint i) const { return ptr[i]; }
    T &operator[](uint i) { return ptr[i]; }

    const T &operator[](const Index2 &id) const
    {
        const uint idx = getLinearIndex(id.x, id.y);
        return ptr[idx];
    }
    T &operator[](const Index2 &id)
    {
        const uint idx = getLinearIndex(id.x, id.y);
        return ptr[idx];
    }

    void clear()
    {
        for (int i = 0, n = size.x * size.y; i < n; i++)
            ptr[i] = T();
    }

private:
    Index2 size;
    T *ptr;
};
};  // namespace asa

#endif
