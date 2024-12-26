#ifndef MNT_POINT_H
#define MNT_POINT_H

#include <fstream> 
#include <proj.h> 

using namespace std;

/* Class used for deserializing a text file containing hydrographic data. */
class Point {

public:
    /** Constructor. */
    Point();
    /** Function allowing conversion from geographic coordinates to Cartesian coordinates.
     * Cartesian coordinates are stored in m_x and m_y.
     * @param P Address of the projection, using Lambert93
     */
    void point_to_pj(PJ* P);
    /** Accessor for Cartesian coordinates.
     * @return double m_x, double m_y
     */
    pair<double, double> coord();
    double m_z; /**< Depth of the point. */

    /** Overload of serialization
     * @param stream
     * @param p
     * @return Displays geographic and Cartesian coordinates if they have been calculated
     */
    friend ostream& operator<<(ostream& stream, const Point& p);
    /** Overload of deserialization
     * @param stream
     * @param p
     * @return Retrieves the latitude, longitude, and depth written on a line
     */
    friend istream& operator>>(istream& stream, Point& p);

private:
    double m_lat, m_long;
    double m_x = 0, m_y = 0;
};


#endif //MNT_POINT_H

