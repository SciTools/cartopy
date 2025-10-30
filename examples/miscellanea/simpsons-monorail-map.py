"""
Recreation of the Monorail Map from The Simpsons
------------------------------------------------

This example demonstrates how to create a minimal outline map of a
defined area of land such as a continent, with optional labels at
specified locations within the region, in the form of a recreation of the
Monorail Map from The Simpsons with humorously oversized labelling (to
imitate handwriting/scribbles) and sparsity of marked locations
(which are all fictional).

Specifically, it aims to recreate to best likeness using Cartopy
the map of pre-Springfield Lyle Lanley Monorail locations from the
iconic episode 'Marge vs. the Monorail' (1993) of the TV Series, as
taken in likeness from the screen grab available at:
https://simpsons.fandom.com/wiki/Brockway.

"""
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader


# Define choices for projection and locations to plot
GEOM_PROJ = ccrs.PlateCarree()
# Not real places, so locations pulled from location on map in still image
# First value 2-tuple is the dot placemap location, second is where to write
# the text label relative to that dot, to best match the map from the show.
LOCATIONS_TO_PLOT = {
    "Ogdenville": [(-111.8, 35.5), (1.5, -2.2), -6],
    "North\nHaverbrook": [(-99.0, 43.5), (2.8, -0.5), -1],
    "Brockway": [(-80.4, 33.6), (-3.4, -1.5), 3],
}


def main():
    # Set up a plot with a thin light blue/white border to make the map inside
    #     Note: for proportions of land mass, font size and border to be as
    #     intended, need to keep 'figsize' and 'dpi' (4:3 ratio) as below.
    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=125, facecolor="white")
    ax.set_facecolor("#AFCBBD")
    ax = fig.add_subplot(111, projection=ccrs.LambertConformal(), frameon=False)
    plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])  # no axes or labels
    ax.set_extent([-120, -72.5, 20, 50], crs=ccrs.Geodetic())  # center on USA

    # Plot only the USA landmass, in a fawn colour with a thin black border
    shpfilename = shpreader.natural_earth(
        resolution="110m", category="cultural", name="admin_0_countries"
    )
    countries = shpreader.Reader(shpfilename).records()
    usa_border = [
        country.geometry
        for country in countries
        if (country.attributes["NAME"] == "United States of America")
    ]
    ax.add_geometries(
        usa_border,
        GEOM_PROJ,
        facecolor="#C39B6A",
        edgecolor="black",
    )

    # From/see example:
    # https://scitools.org.uk/cartopy/docs/latest/gallery/
    # scalar_data/eyja_volcano.html
    # Now add the location labels
    #
    # Use the cartopy interface to create a matplotlib transform object
    # for the Geodetic coordinate system. We will use this along with
    # matplotlib's offset_copy function to define a coordinate system which
    # translates the text by 25 pixels to the left.
    geodetic_transform = GEOM_PROJ._as_mpl_transform(ax)
    text_transform = offset_copy(geodetic_transform, units="dots", x=-25)
    for loc_name, loc_details in LOCATIONS_TO_PLOT.items():
        loc_coords, rel_text_pos, text_rot = loc_details
        ax.plot(
            *loc_coords,
            marker="o",
            color="black",
            markersize=6,
            transform=GEOM_PROJ,
        )

        # Adjust position of location name text relative to location marker
        text_loc_coords = (
            loc_coords[0] + rel_text_pos[0],
            loc_coords[1] + rel_text_pos[1],
        )
        # Text in uppercase, very bold handwriting-like font, as per the
        # screen grab of the map from the show
        ax.text(
            *text_loc_coords,
            loc_name.upper(),
            verticalalignment="center",
            horizontalalignment="left",
            transform=text_transform,
            fontname="Charcoal",  # ensure you have this font available
            fontweight="black",
            fontsize=28,
            rotation=text_rot,  # slightly wonky text for handwritten effect
        )

    leg_text = (
        "Pre-Springfield Lanley\nMonorail locations in TV's\nThe Simpsons\n"
        "(recreation of map at\nsimpsons.fandom.com/\nwiki/Brockway)"
    )

    # Add the 'compass' legend
    ax.text(
        0.14,
        0.10,
        leg_text,
        transform=ax.transAxes,
        fontsize=11,
        horizontalalignment="center",
        verticalalignment="center",
        style="italic",
        bbox=dict(facecolor="#A5B5CE"),
    )

    # Make border symmetrical since default 'rc' file has asymmetric side pad
    fig.tight_layout()
    fig.subplots_adjust(left=0.035, bottom=0.035, right=0.965, top=0.965)
    plt.show()


if __name__ == "__main__":
    main()
