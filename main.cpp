
#pragma comment(lib, "freeglut")
#pragma comment(lib, "glew32")


// SchwarzschildGeodesicExact.cpp
#include <GL/glut.h>
#include <cmath>
#include <vector>
#include <iostream>

const double G = 6.67430e-11; // m^3 kg^-1 s^-2
const double c = 299792458.0; // m/s
const double M = 1.98847e30;

// State variables
double r = 6.98169e10;        // [m]
double pr = 0.0;             // radial velocity at apoapsis [m/s]
double phi = 0.0;            // initial angle [rad]
double pphi = 38860.0 / r;    // initial angular velocity [rad/s] (dφ/dτ), set manually

// Conserved quantities
double E;                  // energy per unit mass (dimensionless)
double L;                  // angular momentum per unit mass [m]

// Visualization
std::vector<float> orbitPoints;
float scale = 1.0f / (1.2f * r); // scaling: Mercury orbit ~70% of window
int windowWidth = 800, windowHeight = 800;

// Schwarzschild radius
const double rs = 2.0 * G * M / (c * c);

// Derivatives
void derivatives(double r, double pr, double& dr_dtau, double& dpr_dtau)
{
    double f = 1.0 - rs / r;

    dr_dtau = pr;

    double dVeff_dr = (G * M) / (r * r) - (L * L) / (r * r * r) + (3.0 * G * M * L * L) / (c * c * r * r * r * r);
    dpr_dtau = -dVeff_dr;
}


void symplecticStep(double dtau)
{
    double f = 1.0 - rs / r;

    // Half-step kick for radial momentum
    double dVeff_dr = (G * M) / (r * r) - (L * L) / (r * r * r) + (3.0 * G * M * L * L) / (c * c * r * r * r * r);
    pr -= 0.5 * dtau * dVeff_dr;

    // Full-step drift for radius
    r += dtau * pr;

    // Update f after radius changes
    f = 1.0 - rs / r;

    // Full-step drift for angle phi (using updated radius)
    phi += (L / (r * r)) * dtau;

    // Recompute potential gradient at new position
    dVeff_dr = (G * M) / (r * r) - (L * L) / (r * r * r) + (3.0 * G * M * L * L) / (c * c * r * r * r * r);

    // Half-step kick for radial momentum
    pr -= 0.5 * dtau * dVeff_dr;

    // Save point for drawing
    float x = static_cast<float>(r * cos(phi) * scale);
    float y = static_cast<float>(r * sin(phi) * scale);
    orbitPoints.push_back(x);
    orbitPoints.push_back(y);
}

void eulerStep(double dtau)
{
    // 1. Compute useful quantity
    double f = 1.0 - rs / r;

    // 2. Compute forces
    double dVeff_dr = (G * M) / (r * r)
        - (L * L) / (r * r * r)
        + (3.0 * G * M * L * L) / (c * c * r * r * r * r);

    // 3. Euler update for momenta
    pr -= dtau * dVeff_dr;

    // 4. Euler update for position
    r += dtau * pr;

    // 5. Euler update for phi (angular motion)
    phi += (L / (r * r)) * dtau;

    // 6. Save point for drawing
    float x = static_cast<float>(r * cos(phi) * scale);
    float y = static_cast<float>(r * sin(phi) * scale);
    orbitPoints.push_back(x);
    orbitPoints.push_back(y);
}


void rk4Step(double dtau)
{
    double k1_r, k1_pr;
    double k2_r, k2_pr;
    double k3_r, k3_pr;
    double k4_r, k4_pr;

    derivatives(r, pr, k1_r, k1_pr);

    derivatives(r + 0.5 * dtau * k1_r, pr + 0.5 * dtau * k1_pr, k2_r, k2_pr);

    derivatives(r + 0.5 * dtau * k2_r, pr + 0.5 * dtau * k2_pr, k3_r, k3_pr);

    derivatives(r + dtau * k3_r, pr + dtau * k3_pr, k4_r, k4_pr);

    r += (dtau / 6.0) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    pr += (dtau / 6.0) * (k1_pr + 2 * k2_pr + 2 * k3_pr + k4_pr);

    // Update phi using conserved angular momentum
    phi += (L / (r * r)) * dtau;

    // Store point
    float x = static_cast<float>(r * cos(phi) * scale);
    float y = static_cast<float>(r * sin(phi) * scale);
    orbitPoints.push_back(x);
    orbitPoints.push_back(y);
}

void initOrbit()
{
    double f = 1.0 - rs / r;

    // Compute conserved quantities exactly from initial velocities
    L = r * r * pphi;

    // From normalization condition, solve for dt/dτ
    double dt_dtau = sqrt(((pr * pr) / f) + (r * r * pphi * pphi) / f + c * c) / c;

    E = f * c * c * dt_dtau;

    std::cout << "Initialized orbit:\n";
    std::cout << "  Schwarzschild radius rs = " << rs << " m\n";
    std::cout << "  E (Energy per mass) = " << E << " J/kg\n";
    std::cout << "  L (Angular momentum per mass) = " << L << " m\n";
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glPointSize(2.0f);

    // Draw black hole
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    glVertex2f(0.0f, 0.0f);
    glEnd();

    // Draw orbit
    glColor3f(1.0, 0.5, 0.0);
    glBegin(GL_LINE_STRIP);
    for (size_t i = 0; i < orbitPoints.size(); i += 2)
        glVertex2f(orbitPoints[i], orbitPoints[i + 1]);

    glEnd();

    glutSwapBuffers();
}

void timer(int value)
{
    //symplecticStep(10000); // simulate proper time
    //rk4Step(10000); // simulate proper time
    eulerStep(10000);

    glutPostRedisplay();
    glutTimerFunc(16, timer, 0);
}

void reshape(int w, int h)
{
    windowWidth = w;
    windowHeight = h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    double aspect = (double)w / (double)h;
    double viewSize = 1.5; // 1.5 units after scaling
    if (aspect >= 1.0)
        glOrtho(-viewSize * aspect, viewSize * aspect, -viewSize, viewSize, -1.0, 1.0);
    else
        glOrtho(-viewSize, viewSize, -viewSize / aspect, viewSize / aspect, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Exact Schwarzschild Geodesic Orbit");

    glClearColor(0.0, 0.0, 0.0, 0.0);

    initOrbit();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, timer, 0);
    glutMainLoop();

    return 0;
}
