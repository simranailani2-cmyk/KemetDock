import py3Dmol

view = py3Dmol.view(width=800, height=600)
view.setStyle({'model': 0}, {})
view.addSurface(py3Dmol.VDW, {'opacity': 0.8, 'color': 'pink'}, {'model': 0})
