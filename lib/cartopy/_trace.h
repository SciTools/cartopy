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


#ifndef _TRACE_H
#define _TRACE_H

#include <iostream>
#include <list>

#include <geodesic.h>
#include <geos_c.h>
#include <proj_api.h>


typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    double x;
    double y;
    double z;
} Vec3;

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
#endif // _TRACE_H
