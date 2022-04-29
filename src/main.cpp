#include "Scene.h"

#include <GL/glut.h>
#include <iostream>

using namespace asa;

// Scene class containing all stuff
Scene gScene;

// Idle function
void idlefunc()
{
    gScene.update();

    glutPostRedisplay();
}

// Keyboard function
void keyboardfunc(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:  // ESC
            exit(0);
            break;
        case 'g':
            gScene.getFluidViz()->toggleVisibleGrid();
            break;
        case 's':
            gScene.pause();
            break;
        case 'p':
            // draw pressure
        default:
            break;
    }
}

// Display function
void displayfunc()
{
    glClear(GL_COLOR_BUFFER_BIT);

    gScene.display();

    glutSwapBuffers();
}

void init()
{
    // Setup the projection matrix.
    glMatrixMode(GL_PROJECTION);
    gluPerspective(45.0, 1.0, 1.0, 10.0);

    // Setup the view and world matrix.
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

// Main program
int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(800, 800);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Basic fluid simulator");

    glutIdleFunc(idlefunc);
    glutKeyboardFunc(keyboardfunc);
    glutDisplayFunc(displayfunc);

    init();
    gScene.init(argc, argv);

    glutMainLoop();

    return 0;
}
