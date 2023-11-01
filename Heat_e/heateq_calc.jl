# import Pkg
using Pkg
# Pkg.add("Configurations")
using Configurations
# Pkg.add("StochasticRounding")
using Plots
using StochasticRounding
using LinearAlgebra
using FileIO
using JLD2

@option struct Parameters
    N           # количество узлов расчетной сетки
    Time        # время интегрирования
    hx          # шаг по пространству
    ht          # шаг по вемени
    TypeNym     # тип чисел с плавающей запятой (Float64, Float32, Float16)
    scheme      # численная схема (RK4, Euler, Euler_KS)
end

p = from_toml(Parameters, "heat_eq.toml")

function RK4(::Type{T},
    w::Vector{Vector{T}},
    tn::Int,
    f,
    h::T,
    DN2) where {T<:AbstractFloat}

    for j = 1:round(Int, tn / h - 1)
        k1 = h * f(DN2, w[j])
        k2 = h * f(DN2, (w[j] + k1) / 2)
        k3 = h * f(DN2, (w[j] + k2) / 2)
        k4 = h * f(DN2, w[j] + k3)

        push!(w, w[j] + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
    end

    return w
end

function Euler(::Type{T}, 
    w::Vector{Vector{T}}, 
    tn::Int, 
    f, 
    h::T,
    DN2) where {T<:AbstractFloat}

    for j = 1:round(Int, tn / h - 1)
        k = h * f(DN2, w[j])

        push!(w, w[j] + k)
    end

    return w
end

function Euler_KS(::Type{T}, 
    w::Vector{Vector{T}}, 
    tn::Int, 
    f, 
    h::T,
    DN2) where {T<:AbstractFloat}
    c = 0
    sum = 0
    for j = 1:round(Int, tn / h - 1)
        y = h * f(DN2, w[j]) .- c
        t = w[j] .+ y
        c = (t .- w[j]) .- y
        sum = t
        push!(w, sum)
    end

    return w
end

function toeplitz(x::Vector{T}) where {T<:Real}
    n = length(x)
    A = zeros(T, n, n)
    for i = 1:n
        for j = 1:n-i+1
            A[i,i+j-1] = x[j]
        end
        for j = n-i+2:n
            A[i, j-(n-i+1)] = x[j]
        end
    end
    return A
 end

function DN_generator(::Type{T}, N, h, Pi) where {T<:AbstractFloat}
    # d = [cot(j * T(h) / 2) / 2 for j = 0:(N - 1)]
    # denom = T.(d)
    h_test = 2 * pi / N
    M = [mod(j, N) != 0 ? (-1) ^ (j+1) * cot(j * h_test / 2) / 2 : 0 for j = 0:(N - 1)]
    DN = transpose.(toeplitz(M))
    return DN
end

function DN2_generator(::Type{T}, N, h, Pi) where {T<:AbstractFloat}
    
    d = [2 * sin(j * h / 2) * sin(j * h / 2) for j = 0:(N - 1)]
    denom = T.(d)
    M = [mod(j, N) != 0 ? (-1) ^ (j+1) / denom[j+1] : (- Pi ^ 2 / (3 * T(h) ^ 2) - 1/6) for j = 0:(N - 1)]
    DN2 = transpose.(toeplitz(M))
    return DN2
end

function DN2_generator_bad(::Type{T}, N, h, Pi) where {T<:AbstractFloat}
    
    d = [2 * sin(j * T(h) / 2) * sin(j * T(h) / 2) for j = 0:(N - 1)]
    # denom = T.(d)
    M = [mod(j, N) != 0 ? (-1) ^ (j+1) / d[j+1] : (- Pi ^ 2 / (3 * T(h) ^ 2) - 1/6) for j = 0:(N - 1)]
    DN2 = transpose.(toeplitz(M))
    return DN2
end

function func(DN2, F)
    return DN2 * F
end

function DN_test(::Type{T}, ht, Pi) where {T<:AbstractFloat}
    M = 100
    errors = zeros(M)
    N = [i*2 for i=1:M]
    for k=1:M
        h = 2 * pi / N[k]
        x = [-pi + i * h for i=1:N[k]]
        f = [exp(sin(2*x[i])) for i = 1:N[k]]
        df = [2*cos(2*x[i]) * exp(sin(2*x[i])) for i = 1:N[k]]
        ddf = [exp(sin(2*x[i]))*(-4*sin(2*x[i])+4*cos(2*x[i])*cos(2*x[i])) for i = 1:N[k]]
        # DN = DN_generator(T, N[k], h, Pi)
        DN = DN2_generator2(N[k], Pi)
        errors[k] = maximum(abs.(DN * f - ddf))
    end

    plot(N, errors, seriestype = :scatter, 
    xaxis=:log, yaxis=:log, 
    title = "Convergence of spectral differentiation", 
    legend = false)
    xlabel!("number of points")
    ylabel!("error")
    savefig("сonvergence_spectral_diffDN2.png")
end

function Heat_equation(::Type{T} = Float64;
    N::Int = p.N,
    Time::Int = p.Time,
    hx::Float64 = p.hx,
    ht::Float64 = p.ht,
    scheme::String = p.sheme) where {T<:AbstractFloat}

    Pi = T.(pi)
    # DN = DN_generator(T, N, hx, Pi)
    # DN2 = DN * DN
    DN2 = DN2_generator_bad(T, N, hx, Pi)

    DN_test(T, ht, Pi)

    # DN216_withoutdenom = DN2_generator2(N, T(hx), Pi)
    # DN264_withoutdenom = DN2_generator2(N, hx, Pi)

    # err1 = findmax(abs.(DN264_withoutdenom .- DN216_withoutdenom))
    # Err_abs1 = err1[1]
    # max1 = abs(DN264_withoutdenom[err1[2]])
    # println("Максимальная ошибка DN2_64-DN216 = ", Err_abs1)

    # err2  = findmax(abs.(DN264_withoutdenom .- DN2))
    # Err_abs2 = err2[1]
    # max2 = abs(DN264_withoutdenom[err2[2]])
    # println("Максимальная ошибка DN2_64-DN2 = ", Err_abs2)

    # Err_otn1 = Err_abs1/max1 *100

    # println("Относительная ошибка DN2_64-DN216 = ", Err_otn1)

    # Err_otn2 = Err_abs2/max2 *100

    # println("Относительная ошибка DN2_16-DN216 = ", Err_otn2)

    x = zeros(T, N)
    x = [T(hx) * j for j=1:N]
    F0 = [cos(x[i]) for i=1:N] # init cond
    ht = ((2 * pi / N) ^ 2) / 5
    # ht = T(ht)
    t = [ht * i for i=1:round(Int, Time/ht)]
    # t = sort(t)
    println(length(t))
    println(length(x))

    if scheme == "RK4"
        return(RK4(T, [F0], Time, func, ht, DN2), x, T.(t))
    elseif scheme == "Euler"
        return(Euler(T, [F0], Time, func, ht, DN2), x, t)
    elseif scheme == "Euler_KS"
        return(Euler_KS(T, [F0], Time, func, ht, DN2), x, t)
    end

end


if p.TypeNym == "Float64"
    FN, x, t = Heat_equation(Float64, N=p.N, Time=p.Time, hx=p.hx, ht=p.ht, scheme=p.scheme)
    # maxeigs = eig(Float64, p.N, p.ht, p.hx)
elseif p.TypeNym == "Float32"
    FN, x, t = Heat_equation(Float32, N=p.N, Time=p.Time, hx=p.hx, ht=p.ht, scheme=p.scheme)
    # maxeigs = eig(Float32, p.N, p.ht, p.hx)
elseif p.TypeNym == "Float16"
    FN, x, t = Heat_equation(Float16, N=p.N, Time=p.Time, hx=p.hx, ht=p.ht, scheme=p.scheme)
    # maxeigs = eig(Float16, p.N, p.ht, p.hx)
end

# FN64, x64, t64 = Heat_equation(Float64, N=p.N, Time=p.Time, hx=p.hx, ht=p.ht, scheme=p.scheme)

M = reduce(vcat, transpose.(FN))
println(length(M))

save("binaryFile.jld2", Dict("M" => M, "x" => x, "t" => t))

function convertMatrixto64(M)
    new_Matrix = Array(fill(0.0, (100, 100)))
    for i in 1:length(M)
        new_Matrix[i] = M[i]
    end
    return new_Matrix
end

function eigvecsvals(::Type{T} = Float64;
    N::Int,
    ht::Float64,
    hx::Float64) where {T<:AbstractFloat}

    Pi = T.(pi)
    DN2 = DN2_generator(T, N, hx, Pi)
    maxvalue = 0
    maxvector = []
    M = I + T.(ht) * DN2
    Matrix = convertMatrixto64(M)
    F = eigen(Matrix)
    val = F.values
    vect = F.vectors
    for i in 1:length(val)
        if abs(val[i]) > maxvalue
            maxvalue = abs(val[i])
        end
    end
    for j in 1:length(val)
        if abs(val[j]) == maxvalue
            maxvector = vect[ : , j ]
        end
    end
    return maxvalue, maxvector
end


function diference(FN1, FN2)
    n = length(FN)
    m = length(FN[1])
    summa = 0
    for i in 1:n
        for j in 1:m
            summa = summa + sum(abs(FN1[i][j] - FN2[i][j]))
        end
    end
    return summa
end

# println(FN[1][1])
# println(diference(FN64, FN64))

# F = eigvecsvals(Float16; p.N, p.ht, p.hx)
# println("Максимальное значение собственного числа = ", F[1])
# println("Собственный вектор, соответствующий этому числу = ", F[2])

# # println(F[:, size(F, 2)])
