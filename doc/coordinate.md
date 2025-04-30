Multiple coordinate systems are used:
- Gate's global (world) system
- Allpix's local (sensor) system

For local, see:
https://allpix-squared.docs.cern.ch/docs/05_geometry_detectors/01_geometry/
=> With local system I use pixel's fraction units

To transform coordinates between systems:
- Use localFractional2globalCoordinates / global2localFractionalCoordinates in utils.py
- See coordinate_transform.jpg