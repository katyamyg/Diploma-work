using Pkg
# Pkg.add("Configurations")
using Configurations

@option struct Parameters
    nu = 1
    L = 35
    N = 100
    t0 = 0
    tN = 300
    dt = 0.05
end

p = from_toml(Parameters, "solver_cfg.toml")

sol_ = Solver_(p.nu, p.L, p.N, p.t0, p.tN, p.dt,) # инициализируем решатель
sol_.solver() # решаем
u = sol_.u    # сохраняем решение
sol_.show()   # выводим решение
