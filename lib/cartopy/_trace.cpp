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
