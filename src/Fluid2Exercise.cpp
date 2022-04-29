#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
namespace
{
float bilinearinterpolation(Vector2 pos, float v00, float v01, float v10, float v11)
{
    float intx = pos.x - floor(pos.x);
    float inty = pos.y - floor(pos.y);
    float vx0 = v01 * intx + v00 * (1 - intx);
    float vx1 = v11 * intx + v10 * (1 - intx);
    return vx1 * inty + vx0 * (1 - inty);
}

float bilinearinterpolation(Vector2 pos, Array2<float> &velarray)
{
    int i = floorf(pos.x), j = floorf(pos.y);
    Index2 s = velarray.getSize();

    float v00 = velarray.getValue(clamp(i, 0, s.x-1), clamp(j, 0, s.y-1));
    float v01 = velarray.getValue(clamp(i+1, 0, s.x-1), clamp(j, 0, s.y-1));
    float v10 = velarray.getValue(clamp(i, 0, s.x-1), clamp(j+1, 0, s.y-1));
    float v11 = velarray.getValue(clamp(i+1, 0, s.x-1), clamp(j+1, 0, s.y-1));

    float vx0 = v00 * (i + 1 - pos.x) + v01 * (pos.x - i);
    float vx1 = v10 * (i + 1 - pos.x) + v11 * (pos.x - i);
    return vx0 * (j + 1 - pos.y) + vx1 * (pos.y - j);
}

Vector3 bilinearinterpolation(Vector2 pos, Array2<Vector3> &inkarray)
{
    int i = floorf(pos.x), j = floorf(pos.y);
    Index2 s = inkarray.getSize();

    Vector3 v00 = inkarray.getValue(clamp(i, 0, s.x-1), clamp(j, 0, s.y-1));
    Vector3 v01 = inkarray.getValue(clamp(i + 1, 0, s.x-1), clamp(j, 0, s.y-1));
    Vector3 v10 = inkarray.getValue(clamp(i, 0, s.x-1), clamp(j + 1, 0, s.y-1));
    Vector3 v11 = inkarray.getValue(clamp(i + 1, 0, s.x-1), clamp(j + 1, 0, s.y-1));

    Vector3 vx0 = v00 * (i + 1 - pos.x) + v01 * (pos.x - i);
    Vector3 vx1 = v10 * (i + 1 - pos.x) + v11 * (pos.x - i);
    return vx0 * (j + 1 - pos.y) + vx1 * (pos.y - j);
}

void acumulate(Array2<float> &velarray, Array2<float> &weight, Vector2 ind, float vel)
{
    float i = floorf(ind.x), j = floorf(ind.y);
    float w;
    if (i >= 0 && i < velarray.getSize().x && j >= 0 && j < velarray.getSize().y) {
        w = (i + 1 - ind.x) * (j + 1 - ind.y);
        velarray.setValue(i, j, velarray.getValue(i, j) + vel * w);
        weight.setValue(i, j, weight.getValue(i, j) + w);
    }

    if (i+1 >= 0 && i + 1 < velarray.getSize().x && j >= 0 && j < velarray.getSize().y) {
        w = (ind.x - i) * (j + 1 - ind.y);
        velarray.setValue(i + 1, j, velarray.getValue(i + 1, j) + vel * w);
        weight.setValue(i + 1, j, weight.getValue(i + 1, j) + w);
    }

    if (i >= 0 && i < velarray.getSize().x && j+1 >= 0 && j + 1 < velarray.getSize().y) {
        w = (i + 1 - ind.x) * (ind.y - j);
        velarray.setValue(i, j + 1, velarray.getValue(i, j + 1) + vel * w);
        weight.setValue(i, j + 1, weight.getValue(i, j + 1) + w);
    }

    if (i+1 >= 0 && i + 1 < velarray.getSize().x && j + 1 >= 0 && j + 1 < velarray.getSize().y) {
        w = (ind.x-i) * (ind.y-j);
        velarray.setValue(i + 1, j + 1, velarray.getValue(i + 1, j + 1) + vel * w);
        weight.setValue(i + 1, j + 1, weight.getValue(i + 1, j + 1) + w);
    }
}

void acumulate(Array2<Vector3> &inkarray, Array2<float> &weight, Vector2 ind, Vector3 ink)
{
    int i = floorf(ind.x), j = floorf(ind.y);
    float w;
    if (i >= 0 && i < inkarray.getSize().x && j >= 0 && j < inkarray.getSize().y) {
        w = (i + 1 - ind.x) * (j + 1 - ind.y);
        inkarray.setValue(i, j, inkarray.getValue(i, j) + ink * w);
        weight.setValue(i, j, weight.getValue(i, j) + w);
    }

    if (i + 1 >= 0 && i + 1 < inkarray.getSize().x && j >= 0 && j < inkarray.getSize().y) {
        w = (ind.x - i) * (j + 1 - ind.y);
        inkarray.setValue(i + 1, j, inkarray.getValue(i + 1, j) + ink * w);
        weight.setValue(i + 1, j, weight.getValue(i + 1, j) + w);
    }

    if (i >= 0 && i < inkarray.getSize().x && j+1 >= 0 && j + 1 < inkarray.getSize().y) {
        w = (i + 1 - ind.x) * (ind.y - j);
        inkarray.setValue(i, j + 1, inkarray.getValue(i, j + 1) + ink * w);
        weight.setValue(i, j + 1, weight.getValue(i, j + 1) + w);
    }

    if ( i + 1 >= 0 && i + 1 < inkarray.getSize().x && j + 1 >= 0 && j + 1 < inkarray.getSize().y) {
        w = (ind.x-i) * (ind.y-j);
        inkarray.setValue(i + 1, j + 1, inkarray.getValue(i + 1, j + 1) + ink * w);
        weight.setValue(i + 1, j + 1, weight.getValue(i + 1, j + 1) + w);
    }
}

inline Vector2 clamp(Vector2 ind, Index2 size)
{
    if (ind.x < 0) {
        ind.x = 0;
    } else if (ind.x > size.x - 1) {
        ind.x = size.x - 1;
    }

    if (ind.y < 0) {
        ind.y = 0;
    } else if (ind.y > size.y - 1) {
        ind.y = size.y - 1;
    }

    return ind;
}

inline Vector2 clamp(Vector2 ind, Vector2 min, Vector2 size)
{
    if (ind.x < min.x) {
        ind.x = min.x;
    } else if (ind.x > size.x) {
        ind.x = size.x;
    }

    if (ind.y < min.y) {
        ind.y = min.y;
    } else if (ind.y > size.y) {
        ind.y = size.y;
    }

    return ind;
}

inline float clamp(float x, float a, float b)
{
    if (x < a) {
        x = a;
    } else if (x > b) {
        x = b;
    }

    return x;
}

inline float random()
{
    return (float)rand() / (float)RAND_MAX;
}

}  // namespace

// init particles
void Fluid2::initParticles()
{
    // particle creation
    /*int total = grid.getSize().x * grid.getSize().y * 4;
    particles.resize(total);
    for (uint i = 0; i < total; i++) {
        float aleatoriox = ((float)rand() / RAND_MAX) * (grid.getDomain().maxPosition.x -
                                                         grid.getDomain().minPosition.x - 2 * grid.getCellDx().x) +
                           grid.getDomain().minPosition.x + grid.getCellDx().x;
        float aleatorioy = ((float)rand() / RAND_MAX) * (grid.getDomain().maxPosition.y -
                                                         grid.getDomain().minPosition.y - 2 * grid.getCellDx().y) +
                           grid.getDomain().minPosition.y + grid.getCellDx().y;
        particles.setPosition(i, Vector2(aleatoriox, aleatorioy));
    }*/

    Vector2 delta = grid.getCellDx()*0.5;
    for (int i = 0; i < grid.getSize().x; i++) {
        for (int j = 0; j < grid.getSize().y; j++) {
            Index2 ind = Index2(i, j);
            Vector2 pos = grid.getCellPos(ind);
            Vector2 auxpos;

            auxpos.x = pos.x + delta.x * random();
            auxpos.y = pos.y + delta.y * random();
            particles.addParticle(auxpos);

            auxpos.x = pos.x - delta.x * random();
            auxpos.y = pos.y + delta.y * random();
            particles.addParticle(auxpos);

            auxpos.x = pos.x + delta.x * random();
            auxpos.y = pos.y - delta.y * random();
            particles.addParticle(auxpos);

            auxpos.x = pos.x - delta.x * random();
            auxpos.y = pos.y - delta.y * random();
            particles.addParticle(auxpos);
        }
    }
}

// advection
void Fluid2::fluidAdvection(const float dt)
{
    if (flipEnabled) {
        // move particles with RK2 with grid velocities

        float vx, vy;
        for (uint i = 0; i < particles.getSize(); i++) {
            Vector2 pos = particles.getPosition(i);

            Vector2 indvx = grid.getFaceIndex(pos, 0);
            vx = bilinearinterpolation(indvx, velocityX);

            Vector2 indvy = grid.getFaceIndex(pos, 1);
            vy = bilinearinterpolation(indvy, velocityY);

            Vector2 vel = Vector2(vx, vy);
            Vector2 midpos = pos + 0.5f * dt * vel;

            indvx = grid.getFaceIndex(midpos, 0);
            vx = bilinearinterpolation(indvx, velocityX);

            indvy = grid.getFaceIndex(pos, 1);
            vy = bilinearinterpolation(indvy, velocityY);

            vel = Vector2(vx, vy);
            pos = pos + dt * vel;

            // ensure particle remains inside the domain

            pos = clamp(pos, grid.getDomain().minPosition, grid.getDomain().maxPosition);

            particles.setPosition(i, pos);
        }

        // create ink grid from particles
        // create velocityX grid from particles
        // create velocityY grid from particles

        Array2<float> inkcoeff(inkRGB.getSize());
        Array2<float> velxcoeff(velocityX.getSize());
        Array2<float> velycoeff(velocityY.getSize());
        velocityX.clear();
        velocityY.clear();
        inkRGB.clear();

        for (uint i = 0; i < particles.getSize(); i++) {
            Vector2 pos = particles.getPosition(i);
            float velx = particles.getVelocity(i).x;
            float vely = particles.getVelocity(i).y;
            Vector3 ink = particles.getInk(i);

            Vector2 indvx = grid.getFaceIndex(pos, 0);
            Vector2 indvy = grid.getFaceIndex(pos, 1);
            Vector2 indcell = grid.getCellIndex(pos);

            acumulate(velocityX, velxcoeff, indvx, velx);
            acumulate(velocityY, velycoeff, indvy, vely);
            acumulate(inkRGB, inkcoeff, indcell, ink);
        }

            // save current state velocities

        for (uint i = 0; i < velocityX.getSize().x; i++) {
            for (uint j = 0; j < velocityX.getSize().y; j++) {
                if (velxcoeff.getValue(i, j) != 0)
                    velocityX.setValue(i, j, velocityX.getValue(i, j) / velxcoeff.getValue(i, j));
            }
        }

        for (uint i = 0; i < velocityY.getSize().x; i++) {
            for (uint j = 0; j < velocityY.getSize().y; j++) {
                if (velycoeff.getValue(i, j) != 0)
                    velocityY.setValue(i, j, velocityY.getValue(i, j) / velycoeff.getValue(i, j));
            }
        }

        for (uint i = 0; i < inkRGB.getSize().x; i++) {
            for (uint j = 0; j < inkRGB.getSize().y; j++) {
                if (inkcoeff.getValue(i, j) != 0)
                    inkRGB.setValue(i, j, inkRGB.getValue(i, j) * (1.f / inkcoeff.getValue(i, j)));
            }
        }

        oldVelocityX.copy(velocityX);
        oldVelocityY.copy(velocityY);

    } else {
        Array2<Vector3> ink(inkRGB.getSize());
        Array2<float> velx(velocityX.getSize());
        Array2<float> vely(velocityY.getSize());
        float vx00, vx01, vx10, vx11, vy00, vy01, vy10, vy11, fvx, fvy;
        float indx0, indx1, indy0, indy1;
        Vector3 ink00, ink01, ink10, ink11, fink;

        for (uint i = 0; i < inkRGB.getSize().x; i++) {
            for (uint j = 0; j < inkRGB.getSize().y; j++) {
                Vector2 pos = grid.getCellPos(Index2(i, j));

                indx1 = min(i + 1, velocityX.getSize().x - 1);
                indy1 = min(j + 1, velocityY.getSize().y - 1);
                Vector2 vel = Vector2((velocityX.getValue(i, j) + velocityX.getValue(indx1, j)) * 0.5f,
                                      (velocityY.getValue(i, j) + velocityY.getValue(i, indy1)) * 0.5f);
                Vector2 auxpos = pos - dt * vel;

                Vector2 indink = grid.getCellIndex(auxpos);

                fink = bilinearinterpolation(indink, inkRGB);
                ink.setValue(i, j, fink);
            }
        }

        for (uint i = 0; i < velocityX.getSize().x; i++) {
            for (uint j = 0; j < velocityX.getSize().y; j++) {
                Vector2 pos = grid.getFaceXPos(Index2(i, j));

                indy0 = max(i - 1, 0u);
                indy1 = min(j + 1, velocityY.getSize().y - 1);
                Vector2 vel = Vector2(velocityX.getValue(i, j),
                                      .25f * (velocityY.getValue(i, j) + velocityY.getValue(indy0, j) +
                                              velocityY.getValue(indy0, indy1) + velocityY.getValue(i, indy1)));
                Vector2 auxpos = pos - dt * vel;

                Vector2 indvx = grid.getFaceIndex(auxpos, 0);

                fvx = bilinearinterpolation(indvx, velocityX);
                velx.setValue(i, j, fvx);
            }
        }

        for (uint i = 0; i < velocityY.getSize().x; i++) {
            for (uint j = 0; j < velocityY.getSize().y; j++) {
                Vector2 pos = grid.getFaceYPos(Index2(i, j));

                indx0 = max(j - 1, 0u);
                indx1 = min(i + 1, velocityX.getSize().x - 1);
                Vector2 vel = Vector2(.25f * (velocityX.getValue(i, j) + velocityX.getValue(i, indx0) +
                                              velocityX.getValue(indx1, j) + velocityX.getValue(indx1, indx0)),
                                      velocityY.getValue(i, j));
                Vector2 auxpos = pos - dt * vel;

                Vector2 indvy = grid.getFaceIndex(auxpos, 1);

                fvy = bilinearinterpolation(indvy, velocityY);
                vely.setValue(i, j, fvy);
            }
        }

        velocityX.copy(velx);
        velocityY.copy(vely);
        inkRGB.copy(ink);
    }
}

void Fluid2::fluidEmission()
{
    if (flipEnabled) {
        // Emitters contribution to particles
        for (uint i = 0; i < particles.getSize(); i++) {
            if (particles.getPosition(i).y <grid.getDomain().minPosition.y*0.9f &&
                particles.getPosition(i).x > -0.1 &&
                particles.getPosition(i).x < 0.1) {
                
                particles.setInk(i,Vector3(1, 0.f, 0.f));
                particles.setVelocity(i,Vector2(0.f, 5.f));
            }
            if (particles.getPosition(i).y < grid.getDomain().minPosition.y * 0.9f &&
                particles.getPosition(i).x > -0.3 &&
                particles.getPosition(i).x < -0.15) {
                particles.setInk(i, Vector3(0.f, 1, 0.f));
                particles.setVelocity(i, Vector2(5.f, 4.f));
            }
            if (particles.getPosition(i).y < grid.getDomain().minPosition.y * 0.9f &&
                particles.getPosition(i).x > 0.15 &&
                particles.getPosition(i).x < 0.3) {
                particles.setInk(i, Vector3(0.f, 0.f, 1));
                particles.setVelocity(i, Vector2(-5.f, 4.f));
            }
        }

    } else {
        // Emitters contribution to grid
        
        for (uint i = grid.getSize().x / 2 - 15; i < grid.getSize().x / 2-10; i++) {
            for (uint j = 0; j < 5; j++) {
                velocityY.setValue(i, j, 4.f);
                velocityX.setValue(i, j, 5.f);

                inkRGB.setValue(i,j, Vector3(0.f, 1.f, 0.f));
            }
        }
        for (uint i = grid.getSize().x / 2 - 2.5; i < grid.getSize().x / 2 + 2.5; i++) {
            for (uint j = 0; j < 5; j++) {
                velocityY.setValue(i, j, 5.f);
                inkRGB.setValue(i, j, Vector3(1.f, 0.f, 0.f));
            }
        }
        for (uint i = grid.getSize().x / 2 +10; i < grid.getSize().x / 2 +15; i++) {
            for (uint j = 0; j < 5; j++) {
                velocityY.setValue(i, j, 4.f);
                velocityX.setValue(i, j, -5.f);
                inkRGB.setValue(i, j, Vector3(0.f, 0.f, 1.f));
            }
        }
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    // Gravity term
     float cons = dt * Scene::kGravity;
    for (uint i = 0; i < velocityY.getSize().x; i++) {
        for (uint j = 0; j < velocityY.getSize().y; j++) {
            velocityY.setValue(i, j, velocityY.getValue(i, j) + cons);
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    // Viscosity term
    float invdx2 = 1 / grid.getCellDx().x / grid.getCellDx().x;
    float invdy2 = 1 / grid.getCellDx().y / grid.getCellDx().y;
    float cons = Scene::kViscosity / Scene::kDensity * dt / invdx2;
    Array2<float> vx(velocityX.getSize());
    Array2<float> vy(velocityY.getSize());
    float aux;

    int nxx = velocityX.getSize().x, nxy = velocityX.getSize().y, nyx = velocityY.getSize().x,
        nyy = velocityY.getSize().y;
    std::vector<double> bx(nxx * nxy);
    std::vector<double> by(nyx * nyy);
    for (int i = 0; i < nxx * nxy; i++) {
        bx[i] = -cons * velocityX[i];
        by[i] = -cons * velocityY[i];
    }
    std::vector<double> fvelx(nxx * nxy);
    std::vector<double> fvely(nyx * nyy);
    multiply(grid.Lx, bx, fvelx);
    multiply(grid.Ly, by, fvely);

    for (int i = 0; i < nxx * nxy; i++) {
        velocityX[i] += fvelx[i];
        velocityY[i] += fvely[i];
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    // Incompressibility / Pressure term

    const Index2 np = pressure.getSize();
    const Index2 sizeP = pressure.getSize();
    const Index2 nx = velocityX.getSize();
    const Index2 ny = velocityY.getSize();

    const Vector2 dx = grid.getCellDx();
    const Vector2 invDx = Vector2(1.0f / dx.x, 1.0f / dx.y);
    const Vector2 invDxPow = Vector2(1.0f / pow(dx.x, 2), 1.0f / pow(dx.y, 2));

    for (int i = 0; i < nx.y; i++) {
        velocityX.setValue(0, i, 0.f);
        velocityX.setValue(nx.x-1, i, 0.f);
    }

    for (int i = 0; i < ny.x; i++) {
        velocityY.setValue(i, 0, 0.f);
        velocityY.setValue(i, ny.y-1, 0.f);
    }

    float pDt = Scene::kDensity / dt;
    std::vector<double> b(np.x * np.y);
    for (unsigned int i = 0; i < np.x; i++) {
        for (unsigned int j = 0; j < np.y; j++) {
            b[pressure.getLinearIndex(i, j)] =
                -pDt * ((velocityX.getValue(i + 1, j) - velocityX.getValue(i, j)) * invDx.x +
                        (velocityY.getValue(i, j + 1) - velocityY.getValue(i, j)) * invDx.y);
        }
    }

    PCGSolver<double> solver;

    double residual_out;
    int iterations_out;
    std::vector<double> p(np.x * np.y);

    solver.set_solver_parameters(1e-6, 1000);
    solver.solve(grid.L, b, p, residual_out, iterations_out);

    for (unsigned int i = 0; i < np.x ; i++) {
        for (unsigned int j = 0; j < np.y; j++) {
            pressure.setValue(i,j, (float)p[pressure.getLinearIndex(i,j)]);
        }
    }

    float dtp = dt / Scene::kDensity;
    for (unsigned int i = 1; i < nx.x - 1; i++) {
        for (unsigned int j = 0; j < nx.y; j++) {
            velocityX.setValue(i, j,
                velocityX.getValue(i, j) - dtp * invDx.x * (pressure.getValue(i, j) - pressure.getValue(i - 1, j)));
        }
    }
    for (unsigned int i = 0; i < ny.x; i++) {
        for (unsigned int j = 1; j < ny.y - 1; j++) {
            velocityY.setValue(i, j,
                velocityY.getValue(i, j) - dtp * (pressure.getValue(i, j) - pressure.getValue(i, j-1)) * invDx.y);
        }
    }

    if (flipEnabled) {

        // PIC/FLIP to update particles velocities

        float vx, vy;
        for (uint i = 0; i < particles.getSize(); i++) {
            Vector2 pos = particles.getPosition(i);

            Vector2 indvx = grid.getFaceIndex(pos, 0);
            vx = bilinearinterpolation(indvx, velocityX);
            Vector2 indvy = grid.getFaceIndex(pos, 1);
            vy = bilinearinterpolation(indvy, velocityY);

            Vector2 newvel = Vector2(vx, vy);

            vx = bilinearinterpolation(indvx, oldVelocityX);
            vy = bilinearinterpolation(indvy, oldVelocityY);

            Vector2 delta = newvel - Vector2(vx, vy);
            Vector2 vel = 0.95f * (particles.getVelocity(i) + delta) + 0.05f * newvel;

            particles.setVelocity(i, vel);
        }
    }
}
}  // namespace asa
