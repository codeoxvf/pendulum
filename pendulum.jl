using Plots

# physical constants
const l₁ = 2.
const l₂ = 1.
const g = 9.8
const m₁ = 2.
const m₂ = 1.

# derivative
# TODO: symbolic diff from energy functions?
const f((θ, θ̇, ϕ, ϕ̇)) = begin
    # coupled Euler-Lagrange equations, solve for θ̈ and ϕ̈ by inverting matrix
    A = [(m₁+m₂)l₁ m₂*l₂*cos(ϕ - θ);
        l₁*cos(ϕ - θ) l₂]
    v = [-m₂ * l₂ * ϕ̇^2 * sin(ϕ - θ) - (m₁ + m₂)g * sin(θ),
        -g * sin(ϕ) - l₁ * θ̇^2 * sin(ϕ - θ)]
    x = A \ v
    [θ̇, x[1], ϕ̇, x[2]]
end

# energy
const U((θ, θ̇, ϕ, ϕ̇)) = begin
    V₁ = m₁ * g * l₁ * (1 - cos(θ))
    V₂ = m₂ * g * (l₁ * (1 - cos(θ)) + l₂ * (1 - cos(ϕ)))
    V₁ + V₂
end

const K((θ, θ̇, ϕ, ϕ̇)) = begin
    T₁ = 0.5m₁ * l₁^2 * θ̇^2
    T₂ = 0.5m₂ * (l₁^2 * θ̇^2 + l₂^2 * ϕ̇^2 + 2l₁ * l₂ * θ̇ * ϕ̇ * cos(ϕ - θ))
    T₁ + T₂
end

const T = 5
const Δt = 0.01
const N = round(Int, T / Δt)
const FPS = 15
const every = round(Int, 1 / FPS / Δt)

# state vector: [θ, θ̇, ϕ, ϕ̇]
y = fill(0., (4, N))

# initial conditions
y[1, 1] = π / 6
y[3, 1] = π / 6

for i in 1:N-1
    # RK4
    k₁ = f(y[:, i])
    k₂ = f(y[:, i] + Δt .* k₁)
    k₃ = f(y[:, i] + Δt / 2 .* k₂)
    k₄ = f(y[:, i] + Δt .* k₃)
    y[:, i+1] = y[:, i] + Δt / 6 .* (k₁ + 2k₂ + 2k₃ + k₄)
end

pendulum_anim = @animate for i in 1:N
    p1 = [l₁ * sin(y[1, i]), -l₁ * cos(y[1, i])]
    p2 = p1 .+ [l₂ * sin(y[3, i]), -l₂ * cos(y[3, i])]

    plot([0, p1[1]], [0, p1[2]], lw=3, marker=:circle, label="")
    plot!([p1[1], p2[1]], [p1[2], p2[2]], lw=3, marker=:circle, label="",
        xlims=(-1.1(l₁ + l₂), 1.1(l₁ + l₂)), ylims=(-1.5(l₁ + l₂), 0),
        aspect_ratio=:equal)
end every every

energy_anim = @animate for i in 1:N
    potential = transpose(mapslices(U, y[:, 1:i], dims=1))
    kinetic = transpose(mapslices(K, y[:, 1:i], dims=1))

    t = Δt * (1:i)
    plot(t, potential, xlabel="Time", ylabel="Energy", label="Potential")
    plot!(t, kinetic, label="Kinetic")
    plot!(t, potential + kinetic, label="Total",
        legend=:bottomright, xlimits=(0, T))

    H₀ = U(y[:, 1]) + K(y[:, 1])
    H = U(y[:, i]) + K(y[:, i])
    annotate!([
        ((0.03, 0.3), ("Initial energy: $H₀", 8, :left)),
        ((0.03, 0.2), ("Current energy: $H", 8, :left))])
end every every

gif(pendulum_anim, "pendulum.gif", fps=FPS)
gif(energy_anim, "energy.gif", fps=FPS)
