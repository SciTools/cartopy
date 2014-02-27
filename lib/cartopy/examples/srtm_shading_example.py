import cartopy.crs as ccrs
from cartopy.io import srtm
import matplotlib.pyplot as plt

plt.figure(figsize=(15, 10))
ax = plt.subplot(111, projection=ccrs.PlateCarree())

elev, crs, extent = srtm.srtm_composite(12, 47, 1, 1)
elev_filled = srtm.fill_srtm(elev, 200.0)
shaded = srtm.shade_srtm(elev_filled, 180.0, 30.0)
plt.imshow(shaded, extent=extent, transform=crs,
           cmap='Greys', origin='lower')
plt.title("SRTM Shaded Relief Map")
gl = ax.gridlines(draw_labels=True,)
gl.xlabels_top = False
gl.ylabels_left = False
plt.show()