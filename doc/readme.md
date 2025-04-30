Gate hits allows to reconstruct ideal cones, as if the Compton camera was perfect. 
Pixel hits from Allpix2 allow to reconstruct cones as would be done with a Timepix3 detector.

A cone is defined with its:
- apex: Compton interaction point
- axis: direction of the scattered photon (vector between Compton and photoelectric interaction points)
- opening angle: angle between the scattered and incident photon directions, calculated from the Compton equation

Compton equation:
cos(opening angle) = 1 - (0.511 * E1_MeV) / (source_MeV * (source_MeV - E1_MeV))
Where:
- E1_MeV = the energy of the scattered photon
- source_MeV = the energy of the incident photon.

# Analysing Gate hits

Geant4 steps can be logged in terminal with sim.g4_verbose, sim.g4_verbose_level_tracking = True, 1 (EventIDs are not logged, hence better do that with small number of events)
- ProcessDefinedStep (Gate) is pre-step, ProcName (G4) is post-step
  => thus when new track are generated in sensor, their ProcessDefinedStep is none.
- KineticEnergy (Gate) is pre-step, KinE (G4) is post-step
- StepLength (Gate) / StepLeng (G4) can be used to match hits (Gate) / steps (G4)

# Analysing Allpix2 pixel hits

Timepix3 can measure particle interactions in frame-based or data driven mode. When used in data driven mode, data is 
stored as a list of pixel hits, where hits contain energy and time.

Allpix2 simulation produces PixelHits objects which resemble pixel hits measured with Timepix3 detectors.
They can be stored as:
- text files using TextWriter module (use `include="PixelHit"` to avoid other data)
- ROOT files using ROOTObjectWriter module

## Pixel hit format
When using TextWriter, Allpix2 stores pixel hits per event with format:
PixelHit pixelID_x, pixelID_y, TOT, local_time, global_time, pixelGlobalCoordX, pixelGlobalCoordY, pixelGlobalCoordZ

## Energy
Timepix3 measures the energy deposited in individual pixels via TOT (Time-Over-Threshold).
When a detector is calibrated with a per-pixel energy calibration procedure, TOT can be converted to energy.

## Time (TOA)

Time-of-Arrival (TOA) is measured with Timepix3 with 1.6 ns granularity.

A so-called time-walk correction can be applied to improve precision, since higher energy deposits induce faster pulses on pixel pre-amplifiers. Precision???

Even though Compton, photo-electric and fluorescent events occur almost simultaneously (within few ps?), the time it takes for the charge carriers to drift to the pixel electrode can be long (depending on semiconductor and bias voltage) and it depends on the depth of interaction.

Drift time = distance / drift speed = distance / (mobility * electric field) 
Drift time = (distance * thickness) / (mobility * voltage)

Drift time depends on the electric field profile in the sensor.
CdTe sensor can have ohmic or Schottky-diode contacts. 
In ohmic contacts, the electric field is constant across the sensor.

Mobility in Silicon (wikipedia):
- electrons: ~1000 cm^2/Vs
- holes: ~450 cm^2/Vs

Mobility in CdTe (wikipedia):
- electrons: ~1100 cm^2/Vs
- holes: ~100 cm^2/Vs

Supposing an ohmic sensor and constant mobility:

Drift time for 1mm in CdTe @ 1000V:
e-: ~10 ns 
holes: ~100 ns

Drift speed in CdTe @ 1000V:
e-: ~100 um/ns
holes: ~10 um/ns

Drift distance during 1.6ns:
e-: ~160 um
hole: ~16 um

## Position

With a single layer Timepix3 camera:
- X/Y coordinates can be calculated from cluster shapes
- Depth of interactions (Z) is more difficult to determine

Delta_Z between two interactions can be calculated from their delta_TOA. See, for example:
10.1088/1748-0221/15/01/C01014

## Clustering hits
When using Allpix2, the output is in pixel hits, as when measuring with a Timepix3 detector. Those are different from gate hits (gHits).
To reconstruct photon interactions in Timepix3 sensor from pixel hits, i.e. their 3D position, time and energy, different algorithms can be used.
