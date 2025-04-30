
# TODO: how to visualize volume sources?
-> use mother volumes for box/sphere?
-> Or as in Gate9: /gate/source/source_name/visualize 1000 yellow 1

# TODO: block below triggers 'WARNING Could not check overlap...' => problem?
pixel = sim.add_volume("Box", "pixel")
pixel.mother, pixel.size = sensor.name, [p, p, thickness]
pixel.material = sensor.material
par = RepeatParametrisedVolume(repeated_volume=pixel)
par.linear_repeat, par.translation = [npix, npix, 1], [p, p, 0]
sim.volume_manager.add_volume(par)

# TODO: deexcitation_ignore_cut impacts number of hits, and depends on cuts

# hits.keep_zero_edep = True # TODO compatible with gHits2cones_byEventID ?

# TODO: cones = pixelClusters2cones(pixel_hits)
