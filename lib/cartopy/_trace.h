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


class Interpolator
{
    public:
    virtual ~Interpolator();
    virtual void set_line(const Point &start, const Point &end);
    virtual Point interpolate(double t)=0;
    virtual Point project(const Point &point)=0;

    protected:
    Point m_start, m_end;
};


class CartesianInterpolator : public Interpolator
{
    public:
    CartesianInterpolator(projPJ src_proj, projPJ dest_proj);
    Point interpolate(double t);
    Point project(const Point &point);

    private:
    projPJ m_src_proj, m_dest_proj;
};


class SphericalInterpolator : public Interpolator
{
    public:
    // XXX Move the constructor and members up to the superclass?
    SphericalInterpolator(projPJ src_proj, projPJ dest_proj);
    void set_line(const Point &start, const Point &end);
    Point interpolate(double t);
    Point project(const Point &point);

    private:
    projPJ m_src_proj, m_dest_proj;
    geod_geodesic m_geod;
    geod_geodesicline m_geod_line;
};


GEOSGeometry *_project_line_string(GEOSContextHandle_t handle,
                                   GEOSGeometry *g_line_string,
                                   Interpolator *interpolator,
                                   GEOSGeometry *g_domain, double threshold);

#endif // _TRACE_H
