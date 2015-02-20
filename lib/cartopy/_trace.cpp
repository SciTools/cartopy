/*
# (C) British Crown Copyright 2010 - 2013, Met Office
#
# This file is part of cartopy.
#
# cartopy is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cartopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with cartopy.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <list>
#include <math.h>
#include <vector>

#include "_trace.h"

#ifdef _MSC_VER
#include <float.h>
#define isnan _isnan
#define isfinite _finite
#endif


Interpolator::~Interpolator(){}

void Interpolator::set_line(const Point &start, const Point &end)
{
    m_start = start;
    m_end = end;
}


CartesianInterpolator::CartesianInterpolator(projPJ src_proj, projPJ dest_proj)
{
    m_src_proj = src_proj;
    m_dest_proj = dest_proj;
}

Point CartesianInterpolator::interpolate(double t)
{
    Point xy;
    xy.x = m_start.x + (m_end.x - m_start.x) * t;
    xy.y = m_start.y + (m_end.y - m_start.y) * t;
    return project(xy);
}

Point CartesianInterpolator::project(const Point &src_xy)
{
    Point dest_xy;
    projLP xy;

    xy.u = src_xy.x;
    xy.v = src_xy.y;

    int status = pj_transform(m_src_proj, m_dest_proj, 1, 1, &xy.u, &xy.v, NULL);
    if (status == -14 || status == -20)
    {
        // -14 => "latitude or longitude exceeded limits"
        // -20 => "tolerance condition error"
        xy.u = xy.v = HUGE_VAL;
    }
    else if (status != 0)
    {
        // TODO: Raise a Python exception instead
        std::cerr << "*******************" << std::endl;
        std::cerr << status << std::endl;
        std::cerr << pj_strerrno(status) << std::endl;
        exit(2);
    }

    dest_xy.x = xy.u;
    dest_xy.y = xy.v;
    return dest_xy;
}


SphericalInterpolator::SphericalInterpolator(projPJ src_proj, projPJ dest_proj)
{
    m_src_proj = src_proj;
    m_dest_proj = dest_proj;
}

void SphericalInterpolator::set_line(const Point &start, const Point &end)
{
    m_start = start;
    m_end = end;

    if (start.x != end.x || start.y != end.y)
    {
        double lon, lat;
        double t, x, y;
        Vec3 end3, axis3;

        // Convert lon/lat to unit vectors
        lon = start.x * DEG_TO_RAD;
        lat = start.y * DEG_TO_RAD;
        t = cos(lat);
        m_start3.x = t * sin(lon);
        m_start3.y = sin(lat);
        m_start3.z = t * cos(lon);

        lon = end.x * DEG_TO_RAD;
        lat = end.y * DEG_TO_RAD;
        t = cos(lat);
        end3.x = t * sin(lon);
        end3.y = sin(lat);
        end3.z = t * cos(lon);

        // Determine the rotation axis for the great circle.
        // axis = ||start x end||
        axis3.x = m_start3.y * end3.z - m_start3.z * end3.y;
        axis3.y = m_start3.z * end3.x - m_start3.x * end3.z;
        axis3.z = m_start3.x * end3.y - m_start3.y * end3.x;
        t = sqrt(axis3.x * axis3.x + axis3.y * axis3.y + axis3.z * axis3.z);
        axis3.x /= t;
        axis3.y /= t;
        axis3.z /= t;

        // Figure out the remaining basis vector.
        // perp = axis x start
        m_perp3.x = axis3.y * m_start3.z - axis3.z * m_start3.y;
        m_perp3.y = axis3.z * m_start3.x - axis3.x * m_start3.z;
        m_perp3.z = axis3.x * m_start3.y - axis3.y * m_start3.x;

        // Derive the rotation angle around the rotation axis.
        x = m_start3.x * end3.x + m_start3.y * end3.y + m_start3.z * end3.z;
        y = m_perp3.x * end3.x + m_perp3.y * end3.y + m_perp3.z * end3.z;
        m_angle = atan2(y, x);
    }
    else
    {
        m_angle = 0.0;
    }
}

Point SphericalInterpolator::interpolate(double t)
{
    Point lonlat;

    if (m_angle == 0.0)
    {
        lonlat = m_start;
    }
    else
    {
        double angle;
        double c, s;
        double x, y, z;
        double lon, lat;

        angle = t * m_angle;
        c = cos(angle);
        s = sin(angle);
        x = m_start3.x * c + m_perp3.x * s;
        y = m_start3.y * c + m_perp3.y * s;
        z = m_start3.z * c + m_perp3.z * s;

        lat = asin(y);
        if(isnan(lat))
        {
            lat = y > 0.0 ? 90.0 : -90.0;
        }
        else
        {
            lat = lat * RAD_TO_DEG;
        }
        lon = atan2(x, z) * RAD_TO_DEG;

        lonlat.x = lon;
        lonlat.y = lat;
    }

    return project(lonlat);
}

Point SphericalInterpolator::project(const Point &lonlat)
{
    Point xy;
    projLP dest;

    //std::cerr << "lon/lat: " << lonlat.x << ", " << lonlat.y;

    dest.u = lonlat.x * DEG_TO_RAD;
    dest.v = lonlat.y * DEG_TO_RAD;

    //std::cerr << " => " << dest.u << ", " << dest.v;

    int status = pj_transform(m_src_proj, m_dest_proj, 1, 1, &dest.u, &dest.v, NULL);
    if (status == -14 || status == -20)
    {
        // -14 => "latitude or longitude exceeded limits"
        // -20 => "tolerance condition error"
        dest.u = dest.v = HUGE_VAL;
    }
    else if (status != 0)
    {
        // TODO: Raise a Python exception instead
        std::cerr << "*******************" << std::endl;
        std::cerr << status << std::endl;
        std::cerr << pj_strerrno(status) << std::endl;
        exit(2);
    }

    //std::cerr << " -> " << dest.u << ", " << dest.v;

    xy.x = dest.u;
    xy.y = dest.v;
    //std::cerr << "xy: " << xy.x << ", " << xy.y << std::endl;
    return xy;
}


typedef std::list<Point> Line;

class LineAccumulator
{
    public:
    LineAccumulator();
    void new_line();
    void add_point(const Point &point);
    void add_point_if_empty(const Point &point);
    GEOSGeometry *as_geom(GEOSContextHandle_t handle);

    std::list<Line>::size_type size() const
    {
        return m_lines.size();
    }

    private:
    std::list<Line> m_lines;
};

LineAccumulator::LineAccumulator()
{
    new_line();
}

void LineAccumulator::new_line()
{
    //std::cerr << "NEW LINE" << std::endl;
    Line line;
    m_lines.push_back(line);
}

void LineAccumulator::add_point(const Point &point)
{
    //std::cerr << "ADD POINT: " << point.x << ", " << point.y << std::endl;
    m_lines.back().push_back(point);
}

void LineAccumulator::add_point_if_empty(const Point &point)
{
    //std::cerr << "ADD POINT IF EMPTY " << m_lines.back().size() << std::endl;
    if (m_lines.back().empty())
    {
        add_point(point);
    }
    //std::cerr << "  FROM EMPTY " << std::endl;
}

bool degenerate_line(const Line &value)
{
    return value.size() < 2;
}

bool close(double a, double b)
{
    return fabs(a - b) <= (1e-8 + 1e-5 * fabs(b));
}

GEOSGeometry *LineAccumulator::as_geom(GEOSContextHandle_t handle)
{
    m_lines.remove_if(degenerate_line);

    if(m_lines.size() > 1)
    {
        //std::cerr << "Checking first & last" << std::endl;
        Point first, last;
        first = m_lines.front().front();
        last = m_lines.back().back();
        //std::cerr << "first: " << first.x << ", " << first.y << std::endl;
        //std::cerr << "last: " << last.x << ", " << last.y << std::endl;
        if(close(first.x, last.x) && close(first.y, last.y))
        {
            m_lines.front().pop_front();
            m_lines.back().splice(m_lines.back().end(), m_lines.front());
            m_lines.pop_front();
        }
    }

    std::vector<GEOSGeometry *> geoms;
    std::list<Line>::const_iterator ilines;
    for(ilines = m_lines.begin(); ilines != m_lines.end(); ++ilines)
    {
        std::list<Point>::const_iterator ipoints;
        int i;

        GEOSCoordSequence *coords = GEOSCoordSeq_create_r(handle, (*ilines).size(), 2);
        for(ipoints = (*ilines).begin(), i = 0; ipoints != (*ilines).end(); ++ipoints, ++i)
        {
            GEOSCoordSeq_setX_r(handle, coords, i, ipoints->x);
            GEOSCoordSeq_setY_r(handle, coords, i, ipoints->y);
        }
        geoms.push_back(GEOSGeom_createLineString_r(handle, coords));
    }

    GEOSGeometry *geom;
    if(geoms.empty())
    {
        geom = GEOSGeom_createEmptyCollection_r(handle, GEOS_MULTILINESTRING);
    }
    else
    {
        geom = GEOSGeom_createCollection_r(handle, GEOS_MULTILINESTRING,
                                           &geoms[0], geoms.size());
    }
    return geom;
}

enum State {
    POINT_IN=1,
    POINT_OUT,
    POINT_NAN
};

State get_state(const Point &point, const GEOSPreparedGeometry *gp_domain,
                GEOSContextHandle_t handle)
{
    State state;

    if (isfinite(point.x) && isfinite(point.y))
    {
        // TODO: Avoid create-destroy
        GEOSCoordSequence *coords = GEOSCoordSeq_create_r(handle, 1, 2);
        GEOSCoordSeq_setX_r(handle, coords, 0, point.x);
        GEOSCoordSeq_setY_r(handle, coords, 0, point.y);
        GEOSGeometry *g_point = GEOSGeom_createPoint_r(handle, coords);
        state = GEOSPreparedCovers_r(handle, gp_domain, g_point) ? POINT_IN : POINT_OUT;
        GEOSGeom_destroy_r(handle, g_point);
    }
    else
    {
        state = POINT_NAN;
    }
    return state;
}

/*
 * Return whether the given line segment is suitable as an
 * approximation of the projection of the source line.
 *
 * t_start: Interpolation parameter for the start point.
 * p_start: Projected start point.
 * t_end: Interpolation parameter for the end point.
 * p_start: Projected end point.
 * interpolator: Interpolator for current source line.
 * threshold: Lateral tolerance in target projection coordinates.
 * handle: Thread-local context handle for GEOS.
 * gp_domain: Prepared polygon of target map domain.
 * inside: Whether the start point is within the map domain.
 */
bool straightAndDomain(double t_start, const Point &p_start,
                       double t_end, const Point &p_end,
                       Interpolator *interpolator, double threshold,
                       GEOSContextHandle_t handle,
                       const GEOSPreparedGeometry *gp_domain,
                       bool inside)
{
    // Straight and in-domain (de9im[7] == 'F')
    
    bool valid;

    // This could be optimised out of the loop.
    if (!(isfinite(p_start.x) && isfinite(p_start.y)))
    {
        valid = false;
    }
    else if (!(isfinite(p_end.x) && isfinite(p_end.y)))
    {
        valid = false;
    }
    else
    {
        // TODO: Re-use geometries, instead of create-destroy!

        // Create a LineString for the current end-point.
        GEOSCoordSequence *coords = GEOSCoordSeq_create_r(handle, 2, 2);
        GEOSCoordSeq_setX_r(handle, coords, 0, p_start.x);
        GEOSCoordSeq_setY_r(handle, coords, 0, p_start.y);
        GEOSCoordSeq_setX_r(handle, coords, 1, p_end.x);
        GEOSCoordSeq_setY_r(handle, coords, 1, p_end.y);
        GEOSGeometry *g_segment = GEOSGeom_createLineString_r(handle, coords);

        // Find the projected mid-point
        double t_mid = (t_start + t_end) * 0.5;
        Point p_mid = interpolator->interpolate(t_mid);

        // Make it into a GEOS geometry
        coords = GEOSCoordSeq_create_r(handle, 1, 2);
        GEOSCoordSeq_setX_r(handle, coords, 0, p_mid.x);
        GEOSCoordSeq_setY_r(handle, coords, 0, p_mid.y);
        GEOSGeometry *g_mid = GEOSGeom_createPoint_r(handle, coords);

        double along = GEOSProjectNormalized_r(handle, g_segment, g_mid);
        if(isnan(along))
        {
            valid = true;
        }
        else
        {
            valid = 0.0 < along && along < 1.0;
            if (valid)
            {
                double separation;
                GEOSDistance_r(handle, g_segment, g_mid, &separation);
                if (inside)
                {
                    // Scale the lateral threshold by the distance from
                    // the nearest end. I.e. Near the ends the lateral
                    // threshold is much smaller; it only has its full
                    // value in the middle.
                    valid = separation <= threshold * 2.0 *
                                            (0.5 - fabs(0.5 - along));
                }
                else
                {
                    // Check if the mid-point makes less than ~11 degree
                    // angle with the straight line.
                    // sin(11') => 0.2
                    // To save the square-root we just use the square of
                    // the lengths, hence:
                    // 0.2 ^ 2 => 0.04
                    double hypot_dx = p_mid.x - p_start.x;
                    double hypot_dy = p_mid.y - p_start.y;
                    double hypot = hypot_dx * hypot_dx + hypot_dy * hypot_dy;
                    valid = ((separation * separation) / hypot) < 0.04;
                }
            }
        }

        if (valid)
        {
            if(inside)
                valid = GEOSPreparedCovers_r(handle, gp_domain, g_segment);
            else
                valid = GEOSPreparedDisjoint_r(handle, gp_domain, g_segment);
        }

        GEOSGeom_destroy_r(handle, g_segment);
        GEOSGeom_destroy_r(handle, g_mid);
    }

    return valid;
}

void bisect(double t_start, const Point &p_start, const Point &p_end,
            GEOSContextHandle_t handle,
            const GEOSPreparedGeometry *gp_domain,
            State &state, Interpolator *interpolator, double threshold,
            double &t_min, Point &p_min, double &t_max, Point &p_max)
{
    double t_current;
    Point p_current;

    // Initialise our bisection range to the start and end points.
    t_min = t_start;
    p_min = p_start;
    t_max = 1.0;
    p_max = p_end;

    // Start the search at the end.
    t_current = t_max;
    p_current = p_max;

    // TODO: See if we can convert the 't' threshold into one based on the
    // projected coordinates - e.g. the resulting line length.
    //
    
    while (fabs(t_max - t_min) > 1.0e-6)
    {
#ifdef DEBUG
        std::cerr << "t: " << t_current << std::endl;
#endif
        bool valid;
        if (state == POINT_IN)
        {
            // Straight and entirely-inside-domain
            valid = straightAndDomain(t_start, p_start, t_current, p_current,
                                      interpolator, threshold,
                                      handle, gp_domain, true);
        }
        else if(state == POINT_OUT)
        {
            // Straight and entirely-outside-domain
            valid = straightAndDomain(t_start, p_start, t_current, p_current,
                                      interpolator, threshold,
                                      handle, gp_domain, false);
        }
        else
        {
            valid = (!isfinite(p_current.x)) || (!isfinite(p_current.y));
        }
#ifdef DEBUG
        std::cerr << "   => valid: " << valid << std::endl;
#endif

        if (valid)
        {
            t_min = t_current;
            p_min = p_current;
        }
        else
        {
            t_max = t_current;
            p_max = p_current;
        }

        t_current = (t_min + t_max) * 0.5;
        p_current = interpolator->interpolate(t_current);
    }
}

void _project_segment(GEOSContextHandle_t handle,
                      const GEOSCoordSequence *src_coords,
                      unsigned int src_idx_from, unsigned int src_idx_to,
                      Interpolator *interpolator,
                      const GEOSPreparedGeometry *gp_domain,
                      double threshold,
                      LineAccumulator &lines)
{
    Point p_current, p_min, p_max, p_end;
    double t_current, t_min, t_max;
    State state;

    GEOSCoordSeq_getX_r(handle, src_coords, src_idx_from, &p_current.x);
    GEOSCoordSeq_getY_r(handle, src_coords, src_idx_from, &p_current.y);
    GEOSCoordSeq_getX_r(handle, src_coords, src_idx_to, &p_end.x);
    GEOSCoordSeq_getY_r(handle, src_coords, src_idx_to, &p_end.y);
#ifdef DEBUG
    std::cerr << "Setting line:" << std::endl;
    std::cerr << "   " << p_current.x << ", " << p_current.y << std::endl;
    std::cerr << "   " << p_end.x << ", " << p_end.y << std::endl;
#endif
    interpolator->set_line(p_current, p_end);
    p_current = interpolator->project(p_current);
    p_end = interpolator->project(p_end);
#ifdef DEBUG
    std::cerr << "Projected as:" << std::endl;
    std::cerr << "   " << p_current.x << ", " << p_current.y << std::endl;
    std::cerr << "   " << p_end.x << ", " << p_end.y << std::endl;
#endif

    t_current = 0.0;
    state = get_state(p_current, gp_domain, handle);

    while(t_current < 1.0 && lines.size() < 500)
    {
        //std::cerr << "Bisecting" << std::endl;
#ifdef DEBUG
        std::cerr << "Working from: " << t_current << " (";
        if (state == POINT_IN)
            std::cerr << "IN";
        else if (state == POINT_OUT)
            std::cerr << "OUT";
        else
            std::cerr << "NAN";
        std::cerr << ")" << std::endl;
        std::cerr << "   " << p_current.x << ", " << p_current.y << std::endl;
        std::cerr << "   " << p_end.x << ", " << p_end.y << std::endl;
#endif
        bisect(t_current, p_current, p_end, handle, gp_domain, state,
               interpolator, threshold,
               t_min, p_min, t_max, p_max);
#ifdef DEBUG
        std::cerr << "   => " << t_min << " to " << t_max << std::endl;
        std::cerr << "   => (" << p_min.x << ", " << p_min.y << ") to (" << p_max.x << ", " << p_max.y << ")" << std::endl;
#endif
        if (state == POINT_IN)
        {
            lines.add_point_if_empty(p_current);
            if (t_min != t_current)
            {
                lines.add_point(p_min);
                t_current = t_min;
                p_current = p_min;
            }
            else
            {
                t_current = t_max;
                p_current = p_max;
                state = get_state(p_current, gp_domain, handle);
                if(state == POINT_IN)
                    lines.new_line();
            }
        }
        else if(state == POINT_OUT)
        {
            if (t_min != t_current)
            {
                t_current = t_min;
                p_current = p_min;
            }
            else
            {
                t_current = t_max;
                p_current = p_max;
                state = get_state(p_current, gp_domain, handle);
                if(state == POINT_IN)
                    lines.new_line();
            }
        }
        else
        {
            t_current = t_max;
            p_current = p_max;
            state = get_state(p_current, gp_domain, handle);
            if(state == POINT_IN)
                lines.new_line();
        }
    }
}

GEOSGeometry *_project_line_string(GEOSContextHandle_t handle,
                                   GEOSGeometry *g_line_string,
                                   Interpolator *interpolator,
                                   GEOSGeometry *g_domain, double threshold)
{
    const GEOSCoordSequence *src_coords = GEOSGeom_getCoordSeq_r(handle, g_line_string);
    unsigned int src_size, src_idx;

    
    const GEOSPreparedGeometry *gp_domain = GEOSPrepare_r(handle, g_domain);

    GEOSCoordSeq_getSize_r(handle, src_coords, &src_size); // check exceptions

    LineAccumulator lines;

    for(src_idx = 1; src_idx < src_size; src_idx++)
    {
        _project_segment(handle, src_coords, src_idx - 1, src_idx,
                      interpolator, gp_domain, threshold, lines);
    }

    GEOSPreparedGeom_destroy_r(handle, gp_domain);

    return lines.as_geom(handle);
}
