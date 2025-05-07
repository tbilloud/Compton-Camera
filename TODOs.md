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

# TODO: include advanced GPU-based reconstruction algorithms, e.g. CoReSi?

# TODO: adapt offline analysis block in main_xxx files to multiple sim runs

# TODO: adapt allpix to opengate 10.0.2

With this version they added a prefix to repeated parametrized volumes. I.e. values in
branch HitUniqueVolumeID values are 'pixel_param-0_0_xxxxx' instead of '0_0_xxxxx' (by default the volume name see to be 0_0). This
prevents the Allpix2 module 'DepositionReader' from reading the volume name properly...
Solutions:
- remove prefix using gate parameters?
- rewrite ROOT file without prefix?
- update detector_name_chars in DepositionReader and change detector name in geometry_conf_content?
  => I tried, looks like theres a problem with character `-` in `pixel_parm-0_0`. 