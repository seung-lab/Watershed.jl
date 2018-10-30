using Watershed
using EMIRT
using HDF5
using Test

aff = h5read(joinpath(@__DIR__, "../assets/piriform.aff.h5"), "main")

#seg = atomicseg(aff)
seg, rg = watershed(aff; is_threshold_relative = true)

#h5write(joinpath(Pkg.dir(), "Watershed/assets/piriform.seg.h5", "seg", seg)
#h5write("seg.h5", "seg", seg)

# compare with segmentation

seg0 = EMIRT.IOs.readseg(joinpath(@__DIR__, "../assets/piriform.seg.h5"))

err = EMIRT.Evaluate.evaluate(seg0, seg)

@show err
@test isapprox(err[:re],  zero(Float64); atol=0.00001)
