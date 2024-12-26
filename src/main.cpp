/** \file main.cpp
 * Main file.
 * To run the program, enter ./create_raster <path to file> <number of pixels>
 * The file must be in a ../files/ folder and have a .txt extension.
 * The number of pixels is an integer.
 */

#include <iostream>
#include <string>

#include <proj.h>
#include <vector>
#include <map>

#include "delaunator.h"
#include "point.h"

using namespace std;

/** Main function. The details of its operation are explained further below.
 *
 * @param filename the name of the file containing the hydrographic data, it must be located in a ../files/ folder
 * @param P projection
 * @param dimension an integer representing the number of pixels along the x-axis
 */

void process(const string& filename, PJ* P, int dimension);

/** Structure containing the 3 integer values R G B defining a color. */
struct Trio {int R, G, B;};
/** Structure containing the 6 doubles defining the coordinates of the 3 points of a triangle. */
struct triangle {double ax, ay, bx, by, cx, cy;};

/** Main function, checks if the correct parameters have been entered.
 *
 * @param argc
 * @param argv
 * @return EXIT_SUCCESS
 */

int main(int argc, char* argv[]) {
    // Initialization of the coordinate reference systems:
    PJ* P = proj_create_crs_to_crs(
            PJ_DEFAULT_CTX,
            "+proj=longlat +datum=WGS84",
            "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs",
            NULL);
    // string filename; //= "../files/Guerledan_Feb19_50cm_wgs84.txt";
    if (argc == 1) cout << "No file provided to process";
    else if (argc == 2) cout << "No dimension provided for the desired image";
    else {
        string filename = argv[1];
        int n = stoi(argv[2]);
        process(filename, P, n);
    }

    return EXIT_SUCCESS;
}

/** Function that checks if a point is inside a triangle with given coordinates
 *
 * @param ai, bi double containing the coordinates of the triangle
 * @param px, py double containing the coordinates of the point
 * @return bool: (px, py) is inside the triangle i
 */

bool in_triangle(double ax,
                 double ay,
                 double bx,
                 double by,
                 double cx,
                 double cy,
                 double px,
                 double py) {
    double X, Y, Z;
    X = (ax - px)*(by - py) - (ay - py)*(bx - px);
    Y = (bx - px)*(cy - py) - (by - py)*(cx - px);
    Z = (cx - px)*(ay - py) - (cy - py)*(ax - px);
    return (X >= 0 && Y >= 0 && Z >= 0) || (X <= 0 && Y <= 0 && Z <= 0);
}

// Choice of the side from which the light comes
double altitude = M_PI/4;
double azimuth = 3*M_PI/4;
/** Function calculating the shading coefficient of a triangle.
 * The function uses the normal to the triangle and arbitrary angles defining the position of the Sun.
 *
 * @param shade double containing the shading coefficient
 * @param xi, yi, zi double containing the coordinates of the 3 points of the triangle
 */

void compute_shade(double &shade, triangle t, double az, double bz, double cz){

    // Normal vectors of the pyramid (1 triangle per plane)
    double nx = (t.by - t.ay) * (cz - az) - (bz - az) * (t.cy - t.ay);
    double ny = (bz - az) * (t.cx - t.ax) - (t.bx - t.ax) * (cz - az);
    double nz = (t.bx - t.ax) * (t.cy - t.ay) - (t.by - t.ay) * (t.cx - t.ax);
    // Complicated formula
    shade = (cos(azimuth) * cos(altitude) * nx + sin(azimuth) * cos(altitude) * ny - sin(altitude) * nz)
            / (sqrt(pow(nx,2) + pow(ny, 2) + pow(nz, 2)) * sqrt(pow(cos(azimuth) * cos(altitude), 2)
                                                                + pow(sin(azimuth) * cos(altitude), 2) + pow(-sin(altitude), 2)));
    // Renormalization
    shade = (1 + shade) / 2;
    // Caps
    if (shade < 0) shade=0;
    if (shade > 1) shade=1;
}

void process(const string& filename, PJ* P, int dimension) {
    // Opening the file and then traversing it
    ifstream in_file(filename);
    if (!in_file.is_open()) throw invalid_argument("Error opening.");

    auto* p = new Point;    /* Converts the coordinates and stores the values */
    auto* v = new vector<double>;    /* Vector provided to the Delaunator */
    map<pair<double, double>, double> prof; /* Depth map based on coordinates */
    double z_min = std::numeric_limits<double>::max();   /* minimum depth */
    double z_max = std::numeric_limits<double>::min();    /* maximum depth */

    cout << "Traversing the file." << endl;
    while (!in_file.eof()) {
        in_file >> *p;
        p->point_to_pj(P);
        v->push_back(p->coord().first);
        v->push_back(p->coord().second);
        prof[p->coord()] = abs(p->m_z);
        if (abs(p->m_z)> z_max) z_max = abs(p->m_z);
        else if (abs(p->m_z)< z_min) z_min = abs(p->m_z);
    }
    in_file.close();

    // Delaunay triangulation and mesh creation
    cout << "Delaunay triangulation and mesh creation." << endl;
    delaunator::Delaunator del(*v);
    map<pair<int, int>, vector<triangle>*> grid; /* map of correspondence between integers and triangles */
    double dx = (del.m_max_x - del.m_min_x)/32; /* gap between two grid cells along x */
    double dy = (del.m_max_y - del.m_min_y)/32; /* gap between two grid cells along y */

    // Definition of a vector for each cell
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 32; j++) {
            grid[pair<int, int> {i, j}] = new vector<triangle>;
        }
    }
   // Traversing and sorting the triangles
    for (int n = 0; n < int (del.triangles.size()); n += 3) {
        triangle t{};
        t.ax = del.coords.at(2 * del.triangles[n]);
        t.ay = del.coords.at(2 * del.triangles[n] + 1);
        t.bx = del.coords.at(2 * del.triangles[n + 1]);
        t.by = del.coords.at(2 * del.triangles[n + 1] + 1);
        t.cx = del.coords.at(2 * del.triangles[n + 2]);
        t.cy = del.coords.at(2 * del.triangles[n + 2] + 1);

        // Coordinates of the bounding box
        double min_x = (t.ax < t.bx) ? t.ax : t.bx;
        min_x = (min_x < t.cx) ? min_x : t.cx;
        double max_x = (t.ax < t.bx) ? t.bx : t.ax;
        max_x = (max_x < t.cx) ? t.cx : max_x;
        double min_y = (t.ay < t.by) ? t.ay : t.by;
        min_y = (min_y < t.cy) ? min_y : t.cy;
        double max_y = (t.ay < t.by) ? t.by : t.ay;
        max_y = (max_y < t.cy) ? t.cy : max_y;

        int n1 = (min_x - del.m_min_x) / dx;
        int n2 = (max_x - del.m_min_x) / dx;
        int n3 = (min_y - del.m_min_y) / dy;
        int n4 = (max_y - del.m_min_y) / dy;

        n2 = n2 >= 32 ? 31 : n2;
        n4 = n4 >= 32 ? 31 : n4;

        // Arbitrary surface value
        double surface = abs((t.bx - t.ax) * (t.cy - t.ay) - (t.cx - t.ax) * (t.by - t.ay)) / 2;
        if (surface < 5) {// Add to the corresponding vectors
            for (int i = n1; i <= n2; i++) {
                for (int j = n3; j <= n4; j++) grid[pair<int, int>{i, j}]->push_back(t);
            }
        }
    }

    cout << "Configuration for file writing." << endl;
    // Creation of the color scale
    ifstream haxby("../src/haxby.txt");   /* file containing the 512 RGB values of the Haxby color scale */
    if (!haxby.is_open()) throw invalid_argument("Error opening haxby.txt.");

    vector<Trio> col;   /* contains the RGB values of the Haxby scale */
    while(!haxby.eof()) {
        Trio RGB{};
        haxby >> RGB.R;
        haxby >> RGB.G;
        haxby >> RGB.B;
        col.push_back(RGB);
    }

    // Determine the image dimensions
    int Nx = dimension; /* dimension entered as an argument */
    int Ny = (del.m_max_y - del.m_min_y) * Nx / (del.m_max_x - del.m_min_x);  /* dimension according to y deduced from Nx */

    // Determine the output file name
    auto* im_name = new string; /* ../images/filename.ppm */
    *im_name = filename;
    im_name->resize(filename.size() - 4);  // Remove the .txt
    im_name->replace(3, 5, "images");   // Change the name of the directory "files" to the new one
    im_name->append("_image.ppm");

    // Opening the new file
    cout << "Opening the new file." << endl;
    ofstream out_file(*im_name); /* binary ppm image corresponding to the input data */
    out_file << "P6" << endl << Nx << " " << Ny << endl << "255" << endl;

    // Estimation of the remaining time
    double percent;
    time_t t0;
    time_t t1;
    time(&t0);
    double chrono;

    // Traversing the pixels
    cout << "Traversing through the pixels." << endl;
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            double px = i * (del.m_max_x - del.m_min_x) / Nx + del.m_min_x; /* abscissa of the pixel in real dimension */
            double py = (Ny - j) * (del.m_max_y - del.m_min_y) / Ny + del.m_min_y;/* ordinate of the pixel in real dimension */

            // Dichotomy search
            double gx = del.m_min_x, gy = del.m_min_y, Lx = del.m_max_x - del.m_min_x, Ly = del.m_max_y - del.m_min_y; /* coordinates of the binary tree node */
            for (int etage = 0; etage < 32; etage ++) {
                if (etage % 2 == 0) {   // Search along x
                    if (px > gx + Lx/2) gx += Lx / 2;
                    Lx = Lx/2;
                }
                else {  // Search along y
                    if (py > gy + Ly/2) gy += Ly/2;
                    Ly = Ly/2;
                }
            }
            int i1 = (gx - del.m_min_x)/dx, i2 = (gy - del.m_min_y)/dy; // integers corresponding to the mesh

            // Search within the triangles
            Trio val{0, 0, 0}; /* color of the pixel */
            double shade = 1; /* shadow of the pixel */
            for (auto t: *grid[pair<int, int> {i1, i2}]) {
                if (in_triangle(t.ax, t.ay, t.bx, t.by, t.cx, t.cy, px, py)) {
                    double az = prof[pair<double, double>(t.ax, t.ay)];
                    double bz = prof[pair<double, double>(t.bx, t.by)];
                    double cz = prof[pair<double, double>(t.cx, t.cy)];
                    double average = (az + bz + cz) / 3;    /* average depth of the three vertices of the triangle */
                    compute_shade(shade, t, az, bz, cz);
                    val = col.at(int(512 * (average - z_min) / (z_max - z_min)));
                    break;
                }
            }

            // Writing to the file
            out_file << (char) (shade * val.R) << (char) (shade * val.G) << (char) (shade * val.B);

            // Estimated remaining time
            percent = 100. * ((i+1) + (j*Nx)) / (Nx * Ny);
            time(&t1);
            chrono = (difftime(t1, t0) * (100 - percent) / percent) / 60;
            cout << "\r" << "[" << percent << "%...] " << "Estimation of remaining time : " << chrono << " minutes." << flush;

        }
    }
    cout << endl;
    // Closing the file
    out_file.close();
}
