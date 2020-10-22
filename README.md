# pandeia-disks
JWST ETC simulations for disks

Quick example:
```
import jwst_pancake.scene as pc_scene
import pandisk as pd
targ = [pd.scene_star('a5v', 4.8)]
targ = pd.add_ring(targ, 20, 5, 1, 180, 15, 45, 30)
pc_scene.rotate_scene(targ, -30)
pd.plot_disk_scene(targ)
```

For install instructions, see https://github.com/spacetelescope/pandeia-coronagraphy
