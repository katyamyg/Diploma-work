using Plots
using JLD2
using FileIO

load("binaryFile.jld2")

heatmap(x, t, M, dpi =300, xlabel="x", ylabel="Time", xtickfontsize=13, ytickfontsize=13, xguidefontsize=13, yguidefontsize=13)
