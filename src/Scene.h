#ifndef SCENE_H
#define SCENE_H

#include "Fluid2.h"
#include "FluidVisualizer2.h"
#include <vector>

namespace asa
{
class Scene
{
public:
    // settings
    static bool pauseFlag;

    static uint nCellsX;
    static uint nCellsY;
    static float step;
    static float kDensity;
    static float kGravity;
    static float kViscosity;
    static uint particlesPerCell;

public:
    Scene();
    ~Scene();

    const Fluid2 *getFluid() const { return fluid; }
    Fluid2 *getFluid() { return fluid; }

    const FluidVisualizer2 *getFluidViz() const { return fluidViz; }
    FluidVisualizer2 *getFluidViz() { return fluidViz; }

    // initialization
    void init();
    void init(int argc, char *argv[]);
    void initAnimation();
    void printSettings();

    // Update
    void pause();
    void update();
    void animate();

    // Display
    void display();

private:
    Fluid2 *fluid;
    FluidVisualizer2 *fluidViz;
};
};  // namespace asa

#endif
