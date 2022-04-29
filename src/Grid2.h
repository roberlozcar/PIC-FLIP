#ifndef GRID_2_H
#define GRID_2_H

#include "Containers/Index2.h"
#include "Math/AABox2.h"
#include "Numeric/SparseMatrix.h"
#include <iostream>
#include <fstream>  

namespace asa
{
class Grid2
{
public:
    Grid2()
        : domain(Vector2(-1, -1), Vector2(1, 1))
        , size(200, 200)
    {
    }

    Grid2(const AABox2 &domain, const Index2 &size) { init(domain, size); }

    void init(const AABox2 &domain, const Index2 &size)
    {
        this->domain = domain;
        this->size = size;
        this->L = CalculateminusLaplacianGridcoef(size.x, size.y);
        this->Lx = CalculateminusLaplacianGridcoef(size.x + 1, size.y);
        this->Ly = CalculateminusLaplacianGridcoef(size.x, size.y + 1);
    }

    AABox2 getDomain() const { return domain; }

    Index2 getSize() const { return size; }

    Vector2 getCellDx() const
    {
        const Vector2 extent(domain.maxPosition - domain.minPosition);
        return Vector2(extent.x / float(size.x), extent.y / float(size.y));
    }

    Index2 getSizeFacesX() const { return Index2(size.x + 1, size.y); }

    Index2 getSizeFacesY() const { return Index2(size.x, size.y + 1); }

    Index2 getSizeNodes() const { return Index2(size.x + 1, size.y + 1); }

    Vector2 getNodePos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x), float(id.y)) * getCellDx());
        return pos;
    }

    Vector2 getCellPos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x) + 0.5f, float(id.y) + 0.5f) * getCellDx());
        return pos;
    }

    Vector2 getFaceXPos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x), float(id.y) + 0.5f) * getCellDx());
        return pos;
    }

    Vector2 getFaceYPos(const Index2 &id) const
    {
        const Vector2 pos(domain.minPosition + Vector2(float(id.x) + 0.5f, float(id.y)) * getCellDx());
        return pos;
    }

    Vector2 getCellIndex(const Vector2 &pos) const
    {
        const Vector2 ispos((pos - domain.minPosition) / getCellDx() - 0.5f);
        return ispos;
    }

    Vector2 getFaceIndex(const Vector2 &pos, const uint axis) const
    {
        Vector2 ispos((pos - domain.minPosition) / getCellDx());
        if (axis == 0)
            ispos.y -= 0.5f;
        if (axis == 1)
            ispos.x -= 0.5f;
        return ispos;
    }

    SparseMatrix<double> CalculateminusLaplacianGridcoef(unsigned int nx, unsigned int ny)
    {
        SparseMatrix<double> A(nx * ny, 5);
        Vector2 invDxPow(1.0f / pow(getCellDx().x, 2), 1.0f / pow(getCellDx().y, 2));
        for (unsigned int i = 0; i < nx; i++) {
            for (unsigned int j = 0; j < ny; j++) {
                int id = j * nx + i;
                if (i > 0) {
                    int id1 = j * nx + i - 1;
                    A.add_to_element(id, id, 1. * invDxPow.x);
                    A.add_to_element(id, id1, -1. * invDxPow.x);
                }
                if (i < nx - 1) {
                    int id1 = j * nx + i + 1;
                    A.add_to_element(id, id, 1. * invDxPow.x);
                    A.add_to_element(id, id1, -1. * invDxPow.x);
                }
                if (j > 0) {
                    int id1 = (j-1) * nx + i;
                    A.add_to_element(id, id, 1. * invDxPow.y);
                    A.add_to_element(id, id1, -1. * invDxPow.y);
                }
                if (j < ny - 1) {
                    int id1 = (j + 1) * nx + i;
                    A.add_to_element(id, id1, -1. * invDxPow.y);
                }
                A.add_to_element(id, id, 1. * invDxPow.y);               
            }
        }
        return A;
    }

public:
    SparseMatrix<double> L;
    SparseMatrix<double> Lx;
    SparseMatrix<double> Ly;

private:
    AABox2 domain;
    Index2 size;
};
};  // namespace asa

#endif
