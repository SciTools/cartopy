/*
# (C) British Crown Copyright 2010 - 2018, Met Office
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
# along with cartopy.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <iostream>
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

    double major_axis, eccentricity_squared;
    pj_get_spheroid_defn(src_proj, &major_axis, &eccentricity_squared);
    geod_init(&m_geod, major_axis, 1 - sqrt(1 - eccentricity_squared));
}

void SphericalInterpolator::set_line(const Point &start, const Point &end)
{
    m_start = start;
    m_end = end;

#if PJ_VERSION > 492
    geod_inverseline(&m_geod_line, &m_geod,
                     m_start.y, m_start.x, m_end.y, m_end.x,
                     GEOD_LATITUDE | GEOD_LONGITUDE);
#else
    double azi1;
    m_a13 = geod_geninverse(&m_geod,
                            m_start.y, m_start.x, m_end.y, m_end.x,
                            NULL, &azi1, NULL, NULL, NULL, NULL, NULL);
    geod_lineinit(&m_geod_line, &m_geod, m_start.y, m_start.x, azi1,
                  GEOD_LATITUDE | GEOD_LONGITUDE);
#endif
}

Point SphericalInterpolator::interpolate(double t)
{
    Point lonlat;

#if PJ_VERSION > 492
    geod_genposition(&m_geod_line, GEOD_ARCMODE, m_geod_line.a13 * t,
                     &lonlat.y, &lonlat.x, NULL, NULL, NULL, NULL, NULL, NULL);
#else
    geod_genposition(&m_geod_line, GEOD_ARCMODE, m_a13 * t,
                     &lonlat.y, &lonlat.x, NULL, NULL, NULL, NULL, NULL, NULL);
#endif

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
