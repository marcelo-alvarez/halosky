import matplotlib.pyplot as plt
import halosky as hs
import numpy as np

hp = hs.halopaint.halopaint(nside=512)

cutout = hp.loadcutout("/global/cfs/cdirs/mp107/exgal/dev/websky-tsz.map")

cutout = np.log10(cutout+1e-15)
plt.imshow(cutout,vmin=-7,vmax=-4,interpolation='nearest')
size = plt.gcf().get_size_inches()
dpi = 4096/size[0] # a bit coarser than native resolution for image of 4096x4096
plt.savefig("websky-tsz.png",dpi=dpi)
plt.show()
