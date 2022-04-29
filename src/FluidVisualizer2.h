#ifndef FLUID_VISUALIZER_2_H
#define FLUID_VISUALIZER_2_H

#include "Fluid2.h"

#include <GL/glut.h>
#include <iostream>

namespace asa
{
class FluidVisualizer2
{
public:
    FluidVisualizer2(const Fluid2 &fluid);
    ~FluidVisualizer2();

    bool getIsVisibleGrid() const { return isVisibleGrid; }

    void init();

    void draw();

    void toggleVisibleGrid();

private:
    void drawGrid();
    void drawInkField();
    void drawVelocityField();
    void drawPressureField();

    void drawParticles();

    void drawBbox(const AABox2 &bbox);

private:
    const Fluid2 &fluid;

    bool isVisibleGrid;

    GLuint texture[3];
};
};  // namespace asa

#endif
