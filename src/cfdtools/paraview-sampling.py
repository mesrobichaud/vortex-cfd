import numpy as np
c = FindSource('Contour1')
c.Isosurfaces = list(np.linspace(-100, 200, 7))
Render()