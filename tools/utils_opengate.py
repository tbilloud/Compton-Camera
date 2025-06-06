# Utility function when using opengate

import os
from pathlib import Path
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

def set_fluorescence(sim):
    sim.physics_manager.global_production_cuts.gamma = 1 * um
    sim.physics_manager.global_production_cuts.electron = 100 * um
    sim.physics_manager.em_parameters.update(
        {'True': True, 'pixe': True, 'deexcitation_ignore_cut': False,
         'auger': True, 'auger_cascade': True})
    # TODO: deexcitation_ignore_cut impacts number of hits, and depends on cuts
    
def get_isotope_data(source):
    """
    Search and print Geant4 radioactive decay data for a given isotope

    Args:
        Z (int): Atomic number
        A (int): Mass number
    """
    import os
    from pathlib import Path

    Z, A = source.ion['Z'], source.ion['A']
    messages = []

    g4_data = os.environ.get('G4RADIOACTIVEDATA')
    if not g4_data:
        messages.append("Warning: G4RADIOACTIVEDATA not set, try to initialize simulation first")
        return "\n".join(messages)

    radioactive_data_path = Path(g4_data)
    isotope_file = f"z{Z}.a{A}"
    possible_files = [
        radioactive_data_path / f"{isotope_file}",
        radioactive_data_path / f"{isotope_file}.z",
        radioactive_data_path / f"{isotope_file}.txt"
    ]

    found_file = None
    for file_path in possible_files:
        if file_path.exists():
            found_file = file_path
            break

    if found_file is None:
        messages.append(f"No decay data found for Z={Z}, A={A}")
        return "\n".join(messages)

    messages.append(f"Found decay data in: {found_file}")

    try:
        with open(found_file, 'r') as f:
            content = f.read()
            messages.append(content)
    except Exception as e:
        messages.append(f"Error reading file: {e}")

    return "\n".join(messages)