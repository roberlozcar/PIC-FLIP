#include "FluidVisualizer2.h"

#include <vector>

namespace asa
{
FluidVisualizer2::FluidVisualizer2(const Fluid2 &fluid)
    : fluid(fluid)
    , isVisibleGrid(false)
{
}

FluidVisualizer2::~FluidVisualizer2()
{
    glDeleteTextures(3, &(texture[0]));
}

void FluidVisualizer2::init()
{
    glGenTextures(3, &(texture[0]));
}

void FluidVisualizer2::draw()
{
    // draw ink RGB field
    drawInkField();

    if (isVisibleGrid) {
        // draw grid
        drawGrid();
        // draw velocity field
        drawVelocityField();
        // draw partices
        drawParticles();
    }
}

void FluidVisualizer2::toggleVisibleGrid()
{
    isVisibleGrid = !isVisibleGrid;
}

void FluidVisualizer2::drawBbox(const AABox2 &bbox)
{
    const Vector2 &a = bbox.minPosition;
    const Vector2 &b = bbox.maxPosition;

    glBegin(GL_LINE_LOOP);
    glColor3f(1.0f, 1.0f, 1.0f);
    glVertex2f(a.x, a.y);
    glVertex2f(a.x, b.y);
    glVertex2f(b.x, b.y);
    glVertex2f(b.x, a.y);
    glEnd();
}

void FluidVisualizer2::drawGrid()
{
    const Vector2 dx = fluid.getGrid().getCellDx();
    const AABox2 gridDomain = fluid.getGrid().getDomain();

    const Vector2 cellCorner1(gridDomain.minPosition);
    const Vector2 cellCorner2(gridDomain.minPosition.x, gridDomain.maxPosition.y - dx.y);
    const Vector2 cellCorner4(gridDomain.maxPosition - dx);
    const Vector2 cellCorner3(gridDomain.maxPosition.x - dx.x, gridDomain.minPosition.y);

    drawBbox(AABox2(gridDomain.minPosition, gridDomain.maxPosition));
    drawBbox(AABox2(cellCorner1, cellCorner1 + dx));
    drawBbox(AABox2(cellCorner2, cellCorner2 + dx));
    drawBbox(AABox2(cellCorner3, cellCorner3 + dx));
    drawBbox(AABox2(cellCorner4, cellCorner4 + dx));
}

void FluidVisualizer2::drawInkField()
{
    const AABox2 domain = fluid.grid.getDomain();
    const Index2 size = fluid.inkRGB.getSize();
    const uint width = size.x;
    const uint height = size.y;

    glBindTexture(GL_TEXTURE_2D, texture[0]);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, width, height, 0, GL_RGB, GL_FLOAT, (GLvoid *)fluid.inkRGB.getData());

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glEnable(GL_TEXTURE_2D);

    glBegin(GL_QUADS);
    glColor3f(0.9f, 0.9f, 0.9f);
    glTexCoord2f(0.0f, 0.0f);
    glVertex2f(domain.minPosition.x, domain.minPosition.y);
    glTexCoord2f(0.0f, 1.0f);
    glVertex2f(domain.minPosition.x, domain.maxPosition.y);
    glTexCoord2f(1.0f, 1.0f);
    glVertex2f(domain.maxPosition.x, domain.maxPosition.y);
    glTexCoord2f(1.0f, 0.0f);
    glVertex2f(domain.maxPosition.x, domain.minPosition.y);
    glEnd();

    glDisable(GL_TEXTURE_2D);
}

void FluidVisualizer2::drawVelocityField()
{
    const float VSCALE = 0.01f;
    const Grid2 &grid = fluid.getGrid();
    const Vector2 dx = grid.getCellDx();

    {
        const Array2<float> &u = fluid.getVelocityX();
        const uint numverts = u.getSize().x * u.getSize().y * 2;
        std::vector<Vector2> verts(numverts);
        for (uint i = 0; i < u.getSize().x; ++i)
            for (uint j = 0; j < u.getSize().y; ++j) {
                const Index2 index(i, j);
                const uint idx = u.getLinearIndex(i, j) * 2;
                Vector2 &head = verts[idx];
                Vector2 &tail = verts[idx + 1];
                head = grid.getFaceXPos(index);
                tail = Vector2(head.x + VSCALE * u[index], head.y);
            }

        glColor3f(0.9f, 0.3f, 0.3f);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, &(verts[0]));
        glDrawArrays(GL_LINES, 0, numverts);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    {
        const Array2<float> &v = fluid.getVelocityY();
        const uint numverts = v.getSize().x * v.getSize().y * 2;
        std::vector<Vector2> verts(numverts);
        for (uint i = 0; i < v.getSize().x; ++i)
            for (uint j = 0; j < v.getSize().y; ++j) {
                const Index2 index(i, j);
                const uint idx = v.getLinearIndex(i, j) * 2;
                Vector2 &head = verts[idx];
                Vector2 &tail = verts[idx + 1];
                head = grid.getFaceYPos(index);
                tail = Vector2(head.x, head.y + VSCALE * v[index]);
            }

        glColor3f(0.3f, 0.9f, 0.3f);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, &(verts[0]));
        glDrawArrays(GL_LINES, 0, numverts);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
}

void FluidVisualizer2::drawParticles(void)
{
    const uint numverts = fluid.getParticles().getSize();

    if (numverts > 0) {
        glColor3f(0.6f, 0.2f, 0.2f);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, &(fluid.getParticles().getPosition(0)));
        glDrawArrays(GL_POINTS, 0, numverts);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
}
}  // namespace asa
