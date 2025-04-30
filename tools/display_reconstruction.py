from napari import view_image, run
from pathlib import Path
import numpy as np
from napari_bbox import BoundingBoxLayer
import matplotlib.pyplot as plt

def display_reconstruction(vol, vsize, vpitch, detector=False):
    # plt.imshow(vol[128,:,:])
    # plt.show()
    # return
    viewer = view_image(vol,
                        translate=tuple(-(v * vpitch) // 2 for v in vsize),
                        axis_labels=['y', 'x', 'z'],
                        scale=[vpitch, vpitch, vpitch],
                        colormap='gray_r')
    viewer.axes.visible = True
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'mm'  # set to None to display no unit
    viewer.scale_bar.length = 10  # length, in units, of the scale bar
    viewer.scale_bar.font_size = 20  # default is 10
    viewer.scale_bar.colored = True  # default value is False
    viewer.scale_bar.color = 'red'  # default value is magenta: (1,0,1,1)
    viewer.scale_bar.position = 'bottom_center'  # default is 'bottom_right'

    # Add detector to the viewer
    # TODO: this makes axes invisible sometimes (according to viewing angle/zoom)
    if detector:
        size, position = detector['size'], detector['position']
        bb_layer = BoundingBoxLayer()
        bb_layer.add([
            [position[0] - size[0] / 2, position[1] - size[1] / 2, position[2] - size[2] / 2],
            [position[0] - size[0] / 2, position[1] - size[1] / 2, position[2] + size[2] / 2],
            [position[0] - size[0] / 2, position[1] + size[1] / 2, position[2] - size[2] / 2],
            [position[0] - size[0] / 2, position[1] + size[1] / 2, position[2] + size[2] / 2],
            [position[0] + size[0] / 2, position[1] - size[1] / 2, position[2] - size[2] / 2],
            [position[0] + size[0] / 2, position[1] - size[1] / 2, position[2] + size[2] / 2],
            [position[0] + size[0] / 2, position[1] + size[1] / 2, position[2] - size[2] / 2],
            [position[0] + size[0] / 2, position[1] + size[1] / 2, position[2] + size[2] / 2]
        ])
        viewer.add_layer(bb_layer)

    run()


if __name__ == "__main__":
    vol = np.load(Path('../output') / "reco_fluoTrue_dopplerFalse.npy")
    display_reconstruction(vol, (256, 256, 256))