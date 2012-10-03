import numpy as np
import matplotlib.pyplot as plt

from cartopy.tests.mpl import image_comparison


@image_comparison(baseline_images=['global_map'])
def test_global_map():
    show = plt.show
    plt.show = lambda *args, **kwargs: None
    import cartopy.examples.global_map as c
    c.main()
    plt.show = show