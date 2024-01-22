import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


data = np.ones((5, 5)) #dummy data
plt.imshow(data, cmap= 'gray', norm = LogNorm(vmin = 1, vmax = 150)) #this is how you logscale matplotlib, vmin and vmax behave like the upper and lower histogram bounds in DS9
plt.show()