using Plots

const l = 1.
const g = 9.8
const m = 1.

# derivative
const f((θ, θ̇)) = [θ̇, -g / l * sin(θ)]

# energy
const U((θ, θ̇)) = m * g * l * (1 - cos(θ))
const K((θ, θ̇)) = 0.5m * l^2 * θ̇^2

const T = 3
const Δt = 0.01
const N = round(Int, T / Δt)

# state vector: [θ, θ̇]
y = fill(0., (2, N))
y[1, 1] = π / 6

for i in 1:N-1
    # RK4
    k₁ = f(y[:, i])
    k₂ = f(y[:, i] + Δt .* k₁)
    k₃ = f(y[:, i] + Δt / 2 .* k₂)
    k₄ = f(y[:, i] + Δt .* k₃)
    y[:, i+1] = y[:, i] + Δt / 6 .* (k₁ + 2k₂ + 2k₃ + k₄)
end

animation = @animate for i in 1:N
    pos = [l * sin(y[1, i]), -l * cos(y[1, i])]
    pendulum_plot = plot([0, pos[1]], [0, pos[2]], xlims=(-1.1l, 1.1l), ylims=(-1.5l, 0),
        aspect_ratio=:equal, lw=3, marker=:circle)

    t = Δt * (1:i)
    potential = transpose(mapslices(U, y[:, 1:i], dims=1))
    kinetic = transpose(mapslices(K, y[:, 1:i], dims=1))
    energy_plot = plot(t, potential, xlabel="Time", ylabel="Energy", label="Potential")
    plot!(t, kinetic, label="Kinetic")
    plot!(t, potential + kinetic, label="Total")

    plot(pendulum_plot, energy_plot, layout=(2, 1))
end

gif(animation, fps=1 / Δt)