using Watershed
using EMIRT
using HDF5

aff = h5read(joinpath(Pkg.dir(), "Watershed/assets/piriform.aff.h5"), "main")

#seg = atomicseg(aff)
seg, rg = watershed(aff)

#h5write(joinpath(Pkg.dir(), "Watershed/assets/piriform.seg.h5", "seg", seg)
#h5write("seg.h5", "seg", seg)

# compare with segmentation

seg0 = readseg(joinpath(Pkg.dir(), "Watershed/assets/piriform.seg.h5"))

err = segerror(seg0, seg)

@show err
@assert err[:re] == 0
