#include <iostream>
#include "point.h"

using namespace std;

Point::Point() :  m_z(), m_lat(), m_long(){};

// Conversion
void Point::point_to_pj(PJ *P) {
    PJ_COORD geo_coord, cartesian_coord;

    geo_coord.lpzt.lam = m_long;
    geo_coord.lpzt.phi = m_lat;
    geo_coord.lpzt.z = m_z;
    cartesian_coord = proj_trans(P, PJ_FWD, geo_coord);

    m_x = cartesian_coord.xy.x;
    m_y = cartesian_coord.xy.y;
}

// Accessor
pair<double, double> Point::coord() {
    pair<double, double> c{m_x, m_y};
    return c;
}

// Serialize
ostream& operator<<(ostream& stream, const Point& p) {
    cout << "point:" << p.m_lat << " " << p.m_long << " " << p.m_z << endl;
    cout << "carte:" << p.m_x << " " << p.m_y << " " << p.m_z << endl;
    return stream;
}

// Deserialize
istream& operator>>(istream& stream, Point& p) {
    stream >> p.m_lat;
    stream >> p.m_long;
    stream >> p.m_z;
    return stream;
}


