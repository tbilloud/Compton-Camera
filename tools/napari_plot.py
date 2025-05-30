import napari
from napari import view_image, run
from napari_bbox import BoundingBoxLayer

def plot_stack_napari(z_slice_stack, vsize):
    """
    Displays a 3D stack of slices using Napari for visualization.
    Compared to matplotlib, it allows for scrolling through the slices

    Parameters:
    -----------
    z_slice_stack : numpy.ndarray or cupy.ndarray
        A 3D array representing the stack of slices to be visualized.
        If using CuPy, it will be converted to a NumPy array.
    vsize : tuple of int
        A tuple representing the size of the volume in the format (depth, height, width).

    Notes:
    ------
    - This function uses Napari to display the stack interactively.
    - The stack is translated so that the center of the volume aligns with the origin.
    - The function assumes that the `xp` module (NumPy or CuPy) is already imported and used for array operations.
    """
    vargs = dict(translate=(-vsize[0] // 2, -vsize[1] // 2), axis_labels=["cone", "x", "y"])
    viewer = napari.view_image(z_slice_stack, **vargs)
    viewer.axes.visible = True
    napari.run()


def plot_reconstruction_napari(vol, vsize, vpitch, detector=False):

    # plt.imshow(vol[128,:,:])
    # plt.show()
    # return

    # if napari:
    #     if xp.__name__ == 'cupy': volume = volume.get()
    #     display_reconstruction(volume, vsize, vpitch, det)

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