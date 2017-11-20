"""
Displaying all 60 zones of the UTM projection
---------------------------------------------

This example displays all 60 zones of the Universal Transverse Mercator
projection next to each other in a figure.

First we create a figure with 60 subplots in one row.
Next we set the projection of each axis in the figure to a specific UTM zone.
Then we add coastlines, gridlines and the number of the zone.
Finally we add a supertitle and display the figure.
"""
__tags__ = ['Lines and polygons']
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


def main():
    # Create a list of integers from 1 - 60
    zones = range(1, 61)

    # Create a figure and array of 60 subplot Axes objects
    fig, axarr = plt.subplots(nrows=1,
                              ncols=len(zones),
                              figsize=(18, 6))

    # Loop through each zone in the list
    for zone in zones:

        # Change the Axes object to a GeoAxes object with a specific UTM zone projection
        axarr[zone - 1] = plt.subplot(1, len(zones), zone,
                                      projection=ccrs.UTM(zone=zone, southern_hemisphere=True))

        # Add coastlines, gridlines and zone number for the subplot
        axarr[zone - 1].coastlines(resolution='110m')
        axarr[zone - 1].gridlines()
        axarr[zone - 1].set_title(zone)

    # Add a supertitle for the figure
    plt.suptitle("UTM Projection - Zones")

    # Display the figure
    plt.show()


if __name__ == '__main__':
    main()
