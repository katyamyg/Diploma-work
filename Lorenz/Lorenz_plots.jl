using Plots
using JLD2
using FileIO
using Pkg
 # Pkg.add("SciPy")
# using SciPy
# Pkg.add(url = "https://github.com/eapax/EarthMover.jl")
using EarthMover
using Configurations

@option struct PlotsParam
    tn
    h
end

p = from_toml(PlotsParam, "Lorenz_plot.toml")

load("example.jld2")

# plt = plot3d(
#     x,
#     y,
#     z,
#     marker = 2,
# )

A64 = Vector{Float64}()
A64_2 = Vector{Float64}()
A32 = Vector{Float32}()
A16 = Vector{Float16}()

for i in 1:length(x64)
    push!(A64, x64[i])
    push!(A64, y64[i])
    push!(A64, z64[i])

    push!(A64_2, x64_2[i])
    push!(A64_2, y64_2[i])
    push!(A64_2, z64_2[i])

    push!(A32, x32[i])
    push!(A32, y32[i])
    push!(A32, z32[i])

    push!(A16, x16[i])
    push!(A16, y16[i])
    push!(A16, z16[i])
end

histogram(A16, xtickfontsize=13, ytickfontsize=13, label = "x,y,z")

# A64 = reshape(A64, 3, length(x64))
# A64_2 = reshape(A64, 3, length(x64))
# A32 = reshape(A64, 3, length(x32))
# A16 = reshape(A64, 3, length(x16))
# println(length(A64))

# WD6_64vs64_2 = []
# WD6_64vs64_2 = [SciPy.stats.wasserstein_distance(A64[1:i], A64_2[1:i]) for i in range(start=3, stop=length(A64), step=3000000)]
# WD = solve_monge(A64[1:3], A64_2[1:3])



# interval1_64 = Vector{Float64}()
# interval2_64 = Vector{Float64}()
# interval3_64 = Vector{Float64}()
# interval4_64 = Vector{Float64}()
# interval5_64 = Vector{Float64}()
# interval6_64 = Vector{Float64}()
# interval7_64 = Vector{Float64}()
# interval8_64 = Vector{Float64}()
# interval9_64 = Vector{Float64}()
# interval10_64 = Vector{Float64}()

# interval1_64_2 = Vector{Float64}()
# interval2_64_2 = Vector{Float64}()
# interval3_64_2 = Vector{Float64}()
# interval4_64_2 = Vector{Float64}()
# interval5_64_2 = Vector{Float64}()
# interval6_64_2 = Vector{Float64}()
# interval7_64_2 = Vector{Float64}()
# interval8_64_2 = Vector{Float64}()
# interval9_64_2 = Vector{Float64}()
# interval10_64_2 = Vector{Float64}()

# interval1_32 = Vector{Float32}()
# interval2_32 = Vector{Float32}()
# interval3_32 = Vector{Float32}()
# interval4_32 = Vector{Float32}()
# interval5_32 = Vector{Float32}()
# interval6_32 = Vector{Float32}()
# interval7_32 = Vector{Float32}()
# interval8_32 = Vector{Float32}()
# interval9_32 = Vector{Float32}()
# interval10_32 = Vector{Float32}()

# interval1_16 = Vector{Float16}()
# interval2_16 = Vector{Float16}()
# interval3_16 = Vector{Float16}()
# interval4_16 = Vector{Float16}()
# interval5_16 = Vector{Float16}()
# interval6_16 = Vector{Float16}()
# interval7_16 = Vector{Float16}()
# interval8_16 = Vector{Float16}()
# interval9_16 = Vector{Float16}()
# interval10_16 = Vector{Float16}()

# num_points = (p.tn / p.h) * 3
# pos = Int(num_points/11)

# for j in (pos):(2*pos)
#     push!(interval1_64, A64[j])
#     push!(interval1_64_2, A64_2[j])
#     push!(interval1_32, A32[j])
#     push!(interval1_16, A16[j])
# end

# for j in (pos):(3*pos)
#     push!(interval2_64, A64[j])
#     push!(interval2_64_2, A64_2[j])
#     push!(interval2_32, A32[j])
#     push!(interval2_16, A16[j])
# end

# for j in (pos):(4*pos)
#     push!(interval3_64, A64[j])
#     push!(interval3_64_2, A64_2[j])
#     push!(interval3_32, A32[j])
#     push!(interval3_16, A16[j])
# end

# for j in (pos):(5*pos)
#     push!(interval4_64, A64[j])
#     push!(interval4_64_2, A64_2[j])
#     push!(interval4_32, A32[j])
#     push!(interval4_16, A16[j])
# end

# for j in (pos):(6*pos)
#     push!(interval5_64, A64[j])
#     push!(interval5_64_2, A64_2[j])
#     push!(interval5_32, A32[j])
#     push!(interval5_16, A16[j])
# end

# for j in (pos):(7*pos)
#     push!(interval6_64, A64[j])
#     push!(interval6_64_2, A64_2[j])
#     push!(interval6_32, A32[j])
#     push!(interval6_16, A16[j])
# end

# for j in (pos):(8*pos)
#     push!(interval7_64, A64[j])
#     push!(interval7_64_2, A64_2[j])
#     push!(interval7_32, A32[j])
#     push!(interval7_16, A16[j])
# end

# for j in (pos):(9*pos)
#     push!(interval8_64, A64[j])
#     push!(interval8_64_2, A64_2[j])
#     push!(interval8_32, A32[j])
#     push!(interval8_16, A16[j])
# end

# for j in (pos):(10*pos)
#     push!(interval9_64, A64[j])
#     push!(interval9_64_2, A64_2[j])
#     push!(interval9_32, A32[j])
#     push!(interval9_16, A16[j])
# end

# for j in (pos):(Int(num_points))
#     push!(interval10_64, A64[j])
#     push!(interval10_64_2, A64_2[j])
#     push!(interval10_32, A32[j])
#     push!(interval10_16, A16[j])
# end

# WD1_64vs64_2 = SciPy.stats.wasserstein_distance(interval1_64, interval1_64_2)
# WD2_64vs64_2 = SciPy.stats.wasserstein_distance(interval2_64, interval2_64_2)
# WD3_64vs64_2 = SciPy.stats.wasserstein_distance(interval3_64, interval3_64_2)
# WD4_64vs64_2 = SciPy.stats.wasserstein_distance(interval4_64, interval4_64_2)
# WD5_64vs64_2 = SciPy.stats.wasserstein_distance(interval5_64, interval5_64_2)
# WD6_64vs64_2 = SciPy.stats.wasserstein_distance(interval6_64, interval6_64_2)
# # WD7_64vs64_2 = SciPy.stats.wasserstein_distance(interval7_64, interval7_64_2)
# # WD8_64vs64_2 = SciPy.stats.wasserstein_distance(interval8_64, interval8_64_2)
# # WD9_64vs64_2 = SciPy.stats.wasserstein_distance(interval9_64, interval9_64_2)
# # WD10_64vs64_2 = SciPy.stats.wasserstein_distance(interval10_64, interval10_64_2)

# t_n = [40000, 60000, 80000, 100000, 120000, 140000]
# # t_n = [40000, 60000, 80000, 100000, 120000, 140000, 160000, 180000, 200000, 220000]

# # plt = plot(t_n, [WD1_64vs64_2,WD2_64vs64_2,WD3_64vs64_2,WD4_64vs64_2,WD5_64vs64_2,WD6_64vs64_2,WD7_64vs64_2,WD8_64vs64_2,WD9_64vs64_2,WD10_64vs64_2], xlabel="t_n", ylabel="WD", label="WD_hi_and_hi", legend=:bottomleft, xaxis=:log, yaxis=:log)
# plt = plot(t_n, [WD1_64vs64_2,WD2_64vs64_2,WD3_64vs64_2,WD4_64vs64_2,WD5_64vs64_2,WD6_64vs64_2], xlabel="t_n", ylabel="WD", label="WD_hi_and_hi", legend=:bottomleft, xaxis=:log, yaxis=:log)