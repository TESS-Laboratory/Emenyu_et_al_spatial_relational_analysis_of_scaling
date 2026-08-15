##running renv
renv::init()
# After installing/updating any package, record the new state:
renv::snapshot()

# Check what's out of sync between renv.lock and your library:
renv::status()

# On a new machine, or to reproduce exactly what's in renv.lock:
renv::restore()