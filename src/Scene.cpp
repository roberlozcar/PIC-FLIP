#include "Scene.h"

#include <GL/glut.h>
#include <iostream>

namespace asa
{
bool Scene::pauseFlag = true;
uint Scene::nCellsX = 200;
uint Scene::nCellsY = 200;
float Scene::step = 0.01f;
float Scene::kDensity = 1.0f;
float Scene::kGravity = -1.0f;
float Scene::kViscosity = 0.001f;
uint Scene::particlesPerCell = 4;

Scene::Scene()
    : fluid(nullptr)
    , fluidViz(nullptr)
{
}

Scene::~Scene()
{
    if (fluidViz)
        delete fluidViz;
    if (fluid)
        delete fluid;
}

void Scene::printSettings()
{
    std::cerr << std::endl << "Current Settings:" << std::endl;
    std::cerr << "\t-gridsize " << nCellsX << " " << nCellsY << std::endl;
    std::cerr << "\t-step " << step << std::endl;
}

void Scene::init(int argc, char *argv[])
{
    int arg = 1;
    while (arg < argc) {
        if (!strcmp(argv[arg], "-gridsize")) {
            arg++;
        } else if (!strcmp(argv[arg], "-step")) {
            arg++;
        } else {
            std::cerr << std::endl << "Unrecognized option " << argv[arg] << std::endl;
            std::cerr << "Usage: practica3.exe -[option1] [settings] -[option2] "
                         "[settings] ..."
                      << std::endl;
            std::cerr << "Options:" << std::endl;
            std::cerr << "\t-test [advection|smoke]" << std::endl;
            std::cerr << "\t-gridsize [gridcells x] [gridcells y]" << std::endl;
            std::cerr << "\t-step [step size in secs]" << std::endl;
            break;
        }
    }

    init();
}

void Scene::init()
{
    printSettings();

    const Index2 gridSize(nCellsX, nCellsY);
    const AABox2 gridDomain(-2.0f, -2.0f, 2.0f, 2.0f);
    const Grid2 grid(gridDomain, gridSize);

    fluid = new Fluid2(grid);
    fluid->init();

    fluidViz = new FluidVisualizer2(*fluid);
    fluidViz->init();

    // initialize test
    initAnimation();
}

void Scene::initAnimation()
{
}

void Scene::pause()
{
    pauseFlag = !pauseFlag;
}

void Scene::update()
{
    if (pauseFlag)
        return;

    animate();
}

void Scene::animate()
{
    fluid->advanceStep(step);
}

void Scene::display()
{
    fluidViz->draw();
}
}  // namespace asa
