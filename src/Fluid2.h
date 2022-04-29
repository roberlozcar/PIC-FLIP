#ifndef FLUID_2_H
#define FLUID_2_H

#include "Grid2.h"
#include "Particles2.h"
#include "Containers/Array2.h"
#include "Math/Vector3.h"

namespace asa
{
class Fluid2
{
public:
    Fluid2(const Grid2 &grid)
        : grid(grid)
        , flipEnabled(false)
    {
    }

    void init()
    {
        velocityX.resize(grid.getSizeFacesX());
        velocityY.resize(grid.getSizeFacesY());
        inkRGB.resize(grid.getSize());
        pressure.resize(grid.getSize());

        if (flipEnabled) {
            oldVelocityX.resize(grid.getSizeFacesX());
            oldVelocityY.resize(grid.getSizeFacesY());

            initParticles();
        }
    }

    const Grid2 &getGrid() const { return grid; }

    const Array2<float> &getVelocityX() const { return velocityX; }

    const Array2<float> &getVelocityY() const { return velocityY; }

    const Array2<float> &getPressure() const { return pressure; }

    const Array2<Vector3> &getInk() const { return inkRGB; }

    const Particles2 &getParticles(void) const { return particles; }

    void advanceStep(const float dt)
    {
        fluidEmission();
        fluidAdvection(dt);
        fluidVolumeForces(dt);
        fluidViscosity(dt);
        fluidPressureProjection(dt);
    }

    // init particles
    void initParticles();

    // emission
    void fluidEmission();

    // advection
    void fluidAdvection(const float dt);

    // external forces
    void fluidVolumeForces(const float dt);

    // viscosity
    void fluidViscosity(const float dt);

    // pressure
    void fluidPressureProjection(const float dt);

public:
    Grid2 grid;

    Array2<float> velocityX;
    Array2<float> velocityY;
    Array2<float> pressure;
    Array2<Vector3> inkRGB;

    bool flipEnabled;

    Array2<float> oldVelocityX;
    Array2<float> oldVelocityY;

    Particles2 particles;
};
};  // namespace asa

#endif
