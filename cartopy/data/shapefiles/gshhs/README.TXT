#	$Id: README.TXT,v 1.21 2011/07/01 20:48:23 guru Exp $

Version 2.2.0 July 2011: The area of small (< 0.1 km^2) polygons
got truncated to 0.  This would cause gshhs to consider them
as lines (borders or rivers) instead of polygons.  Furthermore,
the areas were recomputed using the WGS-84 ellipsoid as the previous
area values were based on a spherical calculation.  Thanks to
José Luis García Pallero for pointing this out.  We now store
the area with a magnitude scale tuned to each polygon.  Also, the
greenwich flag is now a 2-bit flag composed of 1 (crosses Greenwich),
2 (crosses Dateline), 3 (both) or 0 (no such crossing).  See gshhs.[ch] for
details.  Finally, the binary gshhs files now store Antarctica in
-180/+180 range so as to avoid a jump when dumped to ASCII.
Also, the WDBII shapefiles  only had the first 3 levels of rivers;
version 2.2.0 has all 11.
--------------------------------------------------------------------

Version 2.1.1 March 2011: Relatively minor fixes to low-resolution
polygons, including editing errors introduced in v 2.1, removing
a few spikes from 4-5 polygons, and fixing Germany-Poland border
near the Baltic Sea.
--------------------------------------------------------------------

Version 2.1 July 2010: Fixes lack of river-lake flag in the binary
and shapefile release.  Shapefile polygons of level = 2 and with a
negative area are river-lakes.  Also include WDBII border and river
data as shapefiles.

--------------------------------------------------------------------
Global Self-consistent Hierarchical High-resolution Shorelines
		version 2.0 July 15, 2009
		
	Distributed under the Gnu Public License
		
This is the README file for the GSHHS Data distribution.  To read
the data you should get the gshhs supplement to GMT, the Generic
Mapping Tools (gmt.soest.hawaii.edu).  GSHHS appear in GMT in a
different, netCDF format optimized for plotting as huge polygons
are not as efficient.  For more information about how the GSHHS
data were processed, see Wessel and Smith, 1996, JGR.
Many thanks to Tom Kratzke, Metron Inc., for patiently testing
many draft versions of GSHHS and reporting inconsistencies such as
erratic data points and crossings.

Version 2.0 differs from the previous version 1.x in the following
ways.

1.  Free from internal and external crossings and erratic spikes
    at all five resolutions.
2.  The original Eurasiafrica polygon has been split into Eurasia
    (polygon # 0) and Africa (polygon # 1) along the Suez canal.
3.  The original Americas polygon has now been split into North
    America (polygon # 2) and South America (polygon # 3) along
    the Panama canal.
4.  Antarctica is now polygon # 4 and Australia is polygon # 5, in
    all the five resolutions.
5.  Fixed numerous problems, including missing islands and lakes
    in the Amazon and Nile deltas.
6.  Flagged "riverlakes" which are the fat part of major rivers so
    they may easily be identified by users.
7.  Determined container ID for all polygons (== -1 for level 1
    polygons) which is the ID of the polygon that contains a smaller
    polygon.
8.  Determined full-resolution ancestor ID for lower res polygons,
    i.e., the ID of the polygon that was reduced to yield the lower-
    res version.
9.  Ensured consistency across resolutions (i.e., a feature that is
    an island at full resolution should not become a lake in low!).
10. Sorted tables on level, then on the area of each feature.
11. Made sure no feature is missing in one resolution but
    present in the next lower resolution.
12. Store both the actual area of the lower-res polygons and the
    area of the full-resolution ancestor so users may exclude fea-
    tures that represent less that a fraction of the original full
    area.

There was some duplication and wrong levels assigned to maritime
political boundaries in the Persian Gulf that has been fixed.

These changes required us to enhance the GSHHS C-structure used to
read and write the data.  As of version 2.0 the header structure is

struct GSHHS {	/* Global Self-consistent Hierarchical High-resolution Shorelines */
	int id;		/* Unique polygon id number, starting at 0 */
	int n;		/* Number of points in this polygon */
	int flag;	/* = level + version << 8 + greenwich << 16 + source << 24 + river << 25 */
	/* flag contains 5 items, as follows:
	 * low byte:	level = flag & 255: Values: 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
	 * 2nd byte:	version = (flag >> 8) & 255: Values: Should be 7 for GSHHS release 7 (i.e., version 2.0)
	 * 3rd byte:	greenwich = (flag >> 16) & 1: Values: Greenwich is 1 if Greenwich is crossed
	 * 4th byte:	source = (flag >> 24) & 1: Values: 0 = CIA WDBII, 1 = WVS
	 * 4th byte:	river = (flag >> 25) & 1: Values: 0 = not set, 1 = river-lake and level = 2
	 */
	int west, east, south, north;	/* min/max extent in micro-degrees */
	int area;	/* Area of polygon in 1/10 km^2 */
	int area_full;	/* Area of original full-resolution polygon in 1/10 km^2 */
	int container;	/* Id of container polygon that encloses this polygon (-1 if none) */
	int ancestor;	/* Id of ancestor polygon in the full resolution set that was the source of this polygon (-1 if none) */
};

Some useful information:

A) To avoid headaches the binary files were written to be big-endian.
   If you use the GMT supplement gshhs it will check for endian-ness and if needed will
   byte swab the data automatically. If not then you will need to deal with this yourself.

B) In addition to GSHHS we also distribute the files with political boundaries and
   river lines.  These derive from the WDBII data set.

C) As to the best of our knowledge, the GSHHS data are geodetic longitude, latitude
   locations on the WGS-84 ellipsoid.  This is certainly true of the WVS data (the coastlines).
   Lakes, riverlakes (and river lines and political borders) came from the WDBII data set
   which may have been on WGS072.  The difference in ellipsoid is way less then the data
   uncertainties.  Offsets have been noted between GSHHS and modern GPS positions.

D) Originally, the gshhs_dp tool was used on the full resolution data to produce the lower
   resolution versions.  However, the Douglas-Peucker algorithm often produce polygons with
   self-intersections as well as create segments that intersect other polygons.  These problems
   have been corrected in the GSHHS lower resolutions over the years.  If you use gshhs_dp to
   generate your own lower-resolution data set you should expect these problems.

E) The shapefiles release was made by formatting the GSHHS data using the extended GMT/GIS
   metadata understood by OGR, then using ogr2ogr to build the shapefiles.  Each resolution
   is stored in its own subdirectory (e.g., f, h, i, l, c) and each level (1-4) appears in
   its own shapefile.  Thus, GSHHS_h_L3.shp contains islands in lakes for the high res
   data. Because of GIS limitations some polygons that straddle the Dateline (including
   Antarctica) have been split into two parts (east and west).

F) The netcdf-formatted coastlines distributed with GMT derives directly from GSHHS; however
   the polygons have been broken into segments within tiles.  These files are not meant
   to be used by users other than via GMT tools (pscoast, grdlandmask, etc).  The latest
   GMT comes with version 2.0.2 of the netcdf files, still based on GSHHS 2.0.

Paul Wessel	Primary contact: pwessel@hawaii.edu
Walter H. F. Smith

Reference:
Wessel, P. and Smith, W.H.F., 1996. A global, self-consistent, hierarchical, high-resolution
   shoreline database. J. Geophys. Res., 101(B4): 8741–8743.
