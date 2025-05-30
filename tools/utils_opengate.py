import numpy as np
from opengate.logger import global_log
from opengate.geometry.volumes import RepeatParametrisedVolume
from opengate.utility import g4_units

um, mm, keV, MeV, deg, Bq, sec = g4_units.um, g4_units.mm, g4_units.keV, g4_units.MeV, g4_units.deg, g4_units.Bq, g4_units.s


def setup_pixels(sim, npix, sensor, pitch, thickness):
    if not sim.visu:  # because 256 x 256 pixels are too heavy for visualization
        pixel = sim.add_volume("Box", "pixel")
        pixel.mother, pixel.size = sensor.name, [pitch, pitch, thickness]
        pixel.material = sensor.material
        par = RepeatParametrisedVolume(repeated_volume=pixel)
        par.linear_repeat, par.translation = [npix, npix, 1], [pitch, pitch, 0]
        sim.volume_manager.add_volume(par)
    else:
        global_log.warning("Simulation was ran without pixels, analysis will not work.")


def theta_phi(sensor, source):
    # TODO: only works if source and sensor have same x,y coordinates
    sensor_position = np.array(sensor.translation)
    source_position = np.array(source.position.translation)
    sensor_size = np.max(sensor.size[0:1])
    distance = np.linalg.norm(sensor_position - source_position) - sensor.size[
        2] / 2
    phi_deg = 180 - np.degrees(np.arctan(sensor_size / (2 * distance)))
    return [phi_deg * deg, 180 * deg], [0, 360 * deg]
