using Pkg
# Pkg.add("FourierAnalysis")
using PyCall
np = pyimport("numpy")
using PyPlot: plt
using Configurations
using FourierTransforms
using SciPy, FourierAnalysis, Plots
# Pkg.add(PackageSpec(url="https://github.com/eapax/EarthMover.jl", rev="master"))
# using EarthMover


@option struct Parameters
    TypeNum = Float64
    TypeCom = ComplexF64
    nu = 1
    L = 35
    N = 100
    t0 = 0
    tN = 1000
    dt = 0.05
end


mutable struct Solver_
    TypeN   # тип числа с плавающей запятой
    TypeC   # тип комплексного числа с плавающей запятой
    Consist # состояние расчета: первичный или вторичный 
    nu      # параметр, имеющий физический смысл в уравнении
    L       # длина отрезка интегрирования по х
    N       # количество узлов расчетной сетки
    t0      # начальный момент времени
    tN      # конечный момент времени
    dt      # шаг по времени
    k       # коэффициенты k для вычисления операторов дифференцирования
    FL      # оператор дифференцирования для линейной части
    FN      # оператор дифференцирования для нелинейной части
    t       # сетка по времени
    x       # сетка по пространству
    u       # решение
    u_hat   # коэффициенты пространства Фурье
    solver  # метод класса, в котором реализован решатель
    show    # метод класса, демонстрирующий решение

    function Solver_(::Type{T}, ::Type{TC}, Cons, nu, L, N, t0, tN, dt) where {T<:AbstractFloat, TC<:Complex}
        this = new()

        this.TypeN = T
        this.TypeC = TC
        this.Consist = Cons
        this.nu = nu
        this.L = L
        this.N = N
        this.t0 = t0
        this.tN = tN
        this.dt = dt

        this.k = np.arange(-N/2, N/2, 1)
        new_k = T.(this.k)
        this.FL = (((2 * T(pi)) / L ) * new_k) .^ 2 - this.nu * (((2 * T(pi)) / L ) * new_k) .^ 4
        this.FN = - T(1 / 2) * ((1im) * ((2 * T(pi)) / L ) * new_k)

        this.solver = function()
            # количество шагов по времени
            nt = Int((this.tN - this.t0) /this.dt) 

            # сетка
            this.t = np.linspace(start=this.t0, stop=this.tN, num=nt)
            this.x = np.linspace(start=0, stop=this.L, num=this.N) 

            # начальные условия
            if (Cons == 1)
                u0 = -np.cos(T.((2 * T(pi) * this.x) / this.L)) + np.sin(T.((4 * T(pi) * this.x) / this.L))
            else 
                u0 = -np.sin(T.((6 * T(pi) * this.x) / this.L)) + np.cos(T.(16 * (T(pi) * this.x) / this.L)) 
            end
            # println("Начальные условия = ", u0)
            # коэффициенты пространства Фурье
            u0_hat = TC.(FourierTransforms.fftshift(FourierTransforms.fft(u0)))
            # отдельно для нелинейной части
            u0_hat2 = TC.(FourierTransforms.fftshift(FourierTransforms.fft(u0.^2)))

            # массивы с решениями   
            u = zeros(T, (this.N, nt))
            u_hat = zeros(TC, (this.N, nt))
            u_hat2 = zeros(TC, (this.N, nt))

            # инициализация массивов
            u[:, 1] = u0
            u_hat[:, 1] = u0_hat
            u_hat2[:, 1] = u0_hat2

            # основной цикл решателя
            for j in range(1,nt-1)
                uhat_curr = u_hat[:, j]
                uhat_curr2 = u_hat2[:,j]
                if j == 1
                    uhat_prev2 = u_hat2[:,1]
                else
                    uhat_prev2 = u_hat2[:,j-1]
                end
                # схема Кранка-Николсона для линейной части и Адамса-Башфорта для нелинейной части
                # таймстепинг в пространстве коэффициентов Фурье
                u_hat[:, j+1] .= ((1 ./ (1 .- (this.dt / 2) * this.FL)) .* ((1 .+ (this.dt / 2) * this.FL) .* uhat_curr + (((3 / 2) * this.FN) .* (uhat_curr2) - ((1 / 2) * this.FN) .* (uhat_prev2)) * this.dt))
                # возврат в физическое пространство
                u[:,j+1] = np.real(np.fft.ifft(np.fft.ifftshift(u_hat[:,j+1])))
                # корректировка коэффициента
                u_hat[:,j+1] = np.fft.fftshift(np.fft.fft(u[:,j+1]))
                # вычисление нелинейной части
                u_hat2[:,j+1] = np.fft.fftshift(np.fft.fft(u[:,j+1].^2))
            end

            this.u = u
            this.u_hat = u_hat
        end

        this.show = function()
            fig, ax = plt.subplots(figsize=(25,6))
            tt, xx = np.meshgrid(this.t, this.x)
            levels = np.arange(-3, 3, 0.01)
            cs = ax.contourf(tt, xx, this.u)
            plt.yticks(fontsize=20)
            plt.xticks(fontsize=20)
            fig.colorbar(cs)
            ax.set_xlabel("t", fontsize=20)
            ax.set_ylabel("x", fontsize=20)
            ax.set_title("Kuramoto-Sivashinsky", fontsize=30)
            fig.savefig("Kuramoto-Sivashinsky.png")
        end

    return this
    end
end

p = from_toml(Parameters, "solver_cfg.toml")

# sol_ = Solver_(p.nu, p.L, p.N, p.t0, p.tN, p.dt, p.TypeNum) # инициализируем решатель

sol_64 = Solver_(Float64, ComplexF64, 1, 1, 35, 100, 0, 20000, 0.5) # инициализируем решатель
sol_64.solver() # решаем
u64 = sol_64.u    # сохраняем решение
u_hat64 = sol_64.u_hat
# sol_64.show()   # выводим решение

sol_64_2 = Solver_(Float64, ComplexF64, 2, 1, 35, 100, 0, 20000, 0.5) # инициализируем решатель
sol_64_2.solver() # решаем
u64_2 = sol_64_2.u    # сохраняем решение
u_hat64_2 = sol_64_2.u_hat
# sol_64.show()   # выводим решение

sol_16 = Solver_(Float16, ComplexF16, 1, 1, 35, 100, 0, 20000, 0.5) # инициализируем решатель
sol_16.solver() # решаем
u16 = sol_16.u    # сохраняем решение
u_hat16 = sol_16.u_hat
sol_16.show()   # выводим решение

# Поиск взвешенного среднего
function spectrum2(u_hat)
    sp = []
    sp = np.mean((np.absolute([u_hat[50:100, :][i, :] * Complex(-1^(i)/100) for i in range(1, 50)])), axis=1)
    return sp
end 


# Взвешенное среднее 
sp64_weigh = spectrum2(u_hat64)
sp64_2_weigh = spectrum2(u_hat64_2)
sp16_weigh = spectrum2(u_hat16)

spectors = []
push!(spectors, sp64_weigh[1:12], sp64_2_weigh[1:12], sp16_weigh[1:12])
spectors = reduce(vcat,transpose.(spectors))
println(spectors[1, :])
# reshape(spectors, length(spectors), 12)
# import scipy.stats as sts

# WD64vs64_2 = [SciPy.stats.wasserstein_distance(sp64_weigh[i], sp64_2_weigh[i]) for i in range(1, 12)]
# WD64vs16 = [SciPy.stats.wasserstein_distance(sp64_weigh[i], sp16_weigh[i]) for i in range(1, 12)]

WD64vs64_2 = [SciPy.stats.wasserstein_distance(sp64_weigh[1:i], sp64_2_weigh[1:i]) for i in range(start=2, stop=12)]
WD64vs16 = [SciPy.stats.wasserstein_distance(sp64_weigh[1:i], sp16_weigh[1:i]) for i in range(start=2, stop=12)]


# t_nn = [i for i=1:10]
# plot(t_nn, sp64_weigh[1:10], label="Float64")
# plot!(t_nn, sp64_2_weigh[1:10], label="Float64_2")
# plot!(t_nn, sp16_weigh[1:10], label="Float16", xlabel="волновое число", ylabel = "средняя амплитуда")

t_nn = [i for i=1:11]
plot(t_nn, WD64vs64_2, label="Float64vsFloat64_2", xlabel="волновое число", ylabel = "расстояние Вассерштейна")
plot!(t_nn, WD64vs16, label="Float64vsFloat16")
# plot!(t_nn, sp16_weigh[1:10], label="Float16", xlabel="индекс частоты", ylabel = "средняя амплитуда")