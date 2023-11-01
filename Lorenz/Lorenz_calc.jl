using Pkg
# Pkg.add("Configurations")
# Pkg.add("Plots")
# Pkg.add("FileIO")
# Pkg.add("JLD2")
using Configurations
using FileIO
using JLD2

@option struct Parameters
    TypeNym
    x0
    tn
    h
    a
    b
    c
    x1
    y1
    z1
    scheme
end

p = from_toml(Parameters, "Lorenz.toml")



# Модель Лоренца, состоящая из трех ОДУ

function f1_16(x, y, a)
    return a * (y - x)
end

function f2_16(x, y, z, b)
    return x * (b - z) - y
end

function f3_16(x, y, z, c)
    return x * y - c * z
end


# Метод интегрирования системы ОДУ - метод Рунге-Кутта 4-го порядка

function RK4(::Type{T},
    x0::Float64,
    tn::Int,
    h::Float64,
    a::Float64,
    b::Float64,
    c::Float64,
    x1::Float64,
    y1::Float64,
    z1::Float64) where {T<:AbstractFloat}

    m = Int64(round((tn - x0) / h))
    x = zeros(T, m + 1)
    y = zeros(T, m + 1)
    z = zeros(T, m + 1)

    x1 = T.(x1)
    y1 = T.(y1)
    z1 = T.(z1)
    a = T.(a)
    b = T.(b)
    c = T.(c)

    x[1] = x1;
    y[1] = y1;
    z[1] = z1;

    for i in 1:m
      kx1 = h * f1_16(x[i], y[i], a)
      ky1 = h * f2_16(x[i], y[i], z[i], b)
      kz1 = h * f3_16(x[i], y[i], z[i], c)

      kx2 = h * f1_16(x[i] + h/2, y[i] + 1/2 * ky1, a)
      ky2 = h * f2_16(x[i] + h / 2, y[i] + 1 / 2 * ky1, z[i] + 1 / 2 * kz1, b)
      kz2 = h * f3_16(x[i] + h / 2, y[i] + 1 / 2 * ky1, z[i] + 1 / 2 * kz1, c)

      kx3 = h * f1_16(x[i] + h/2, y[i] + 1/2 * ky2, a)
      ky3 = h * f2_16(x[i] + h / 2, y[i] + 1 / 2 * ky2, z[i] + 1 / 2 * kz2, b)
      kz3 = h * f3_16(x[i] + h / 2, y[i] + 1 / 2 * ky2, z[i] + 1 / 2 * kz2, c)

      kx4 = h * f1_16(x[i] + h, y[i] + ky3, a)
      ky4 = h * f2_16(x[i] + h, y[i] + ky3, z[i] + kz3, b)
      kz4 = h * f3_16(x[i] + h, y[i] + ky3, z[i] + kz3, c)

      x[i + 1] = x[i] + (1 / 6) * (kx1 + 2 * kx2 + 2 * kx3 + kx4)
      y[i + 1] = y[i] + (1 / 6) * (ky1 + 2 * ky2 + 2 * ky3 + ky4)
      z[i + 1] = z[i] + (1 / 6) * (kz1 + 2 * kz2 + 2 * kz3 + kz4)
    end

      return x, y, z
end


function Lorenz(::Type{T} = Float64;
    x0::Float64 = p.x0,
    tn::Int = p.tn,
    h::Float64 = p.h,
    a::Float64 = p.a,
    b::Float64 = p.b,
    c::Float64 = p.c,
    x1::Float64 = p.x1,
    y1::Float64 = p.y1,
    z1::Float64 = p.z1,
    scheme::String = p.scheme) where {T<:AbstractFloat}

    if scheme == "RK4"
        return(RK4(T, x0, tn, h, a, b, c, x1, y1, z1))
    else 
        throw(error("Other schemes not implemented yet."))
    end
    
end


# if p.TypeNym == "Float64"
#     x, y, z = Lorenz(Float64, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")
# elseif p.TypeNym == "Float32"
#     x, y, z = Lorenz(Float32, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")
# elseif p.TypeNym == "Float16"
#     x, y, z = Lorenz(Float16, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")
# end

x64, y64, z64 = Lorenz(Float64, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")
x64_2, y64_2, z64_2 = Lorenz(Float64, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=0.9, y1=0.3, z1=2.1, scheme="RK4")
x32, y32, z32 = Lorenz(Float32, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")
x16, y16, z16 = Lorenz(Float16, x0=p.x0, tn=p.tn, h=p.h, a=p.a, b=p.b, c=p.c, x1=p.x1, y1=p.y1, z1=p.z1, scheme="RK4")

save("example.jld2", Dict("x64" => x64, "y64" => y64, "z64" => z64, "x64_2" => x64_2, "y64_2" => y64_2, "z64_2" => z64_2, "x32" => x32, "y32" => y32, "z32" => z32, "x16" => x16, "y16" => y16, "z16" => z16))