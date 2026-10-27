using GLMakie, ColorSchemes, ColorTypes, MakieThemes

c = 299792458

# Specifications
lambda0 = 1550e-9       # m
f0 = c / lambda0        # Hz
channel_spacing = 100e9 # Hz

# Ring parameters
R = 10e-6           # Radius (m)
ng = 4.3            # Group index
L = 2 * pi * R   # Ring circumference

# Free spectral range approximation
fsr = c / (ng * L)
fsr_w = lambda0^2 / (ng * L)

println("Target frequency: ", f0 / 1e12, " THz")
println("Ring FSR: ", fsr / 1e12, " THz")
println("Ring FSR wavelength: ", fsr_w * 1e9, " nm")
println("FSR / channel spacing: ", fsr / channel_spacing)

# Coupling and loss parameters
a = 0.99
t1 = 0.975
t2 = 0.97

# Through port transmission
function ring_through(f, R, ng, a, t1, t2, f0)
    L = 2 * pi * R

    # Phase relative to the selected resonance
    phi = 2 * pi * ng * L * (f - f0) / c

    denominator = 1 + (a * t1 * t2)^2 - 2 * a * t1 * t2 * cos(phi)

    T_through = (t1^2 + (a * t2)^2 - 2 * a * t1 * t2 * cos(phi)) / denominator

    return T_through
end

# Drop port transmission
function ring_drop(f, R, ng, a, t1, t2, f0)
    L = 2 * pi * R

    # Phase relative to the selected resonance
    phi = 2 * pi * ng * L * (f - f0) / c

    denominator = 1 + (a * t1 * t2)^2 - 2 * a * t1 * t2 * cos(phi)

    T_drop = a * (1 - t1^2) * (1 - t2^2) / denominator

    return T_drop
end

span = 1.5 * fsr
f_range = LinRange(f0 - span, f0 + span, 20001)

T_through = [ring_through(f, R, ng, a, t1, t2, f0) for f in f_range]
T_drop = [ring_drop(f, R, ng, a, t1, t2, f0) for f in f_range]

fig = Figure(fontsize=36, px_per_unit = 2, size= 3 .*(500, 170), figure_padding = 25)
ax = Axis(fig[1, 1], xlabel = L"f \ (THz)", 
    ylabel = L"T \ (dB)", aspect=4)
lines!(ax, f_range / 1e12, 10*log.(T_through), color="blue", label="through port")
lines!(ax, f_range / 1e12, 10*log.(T_drop), color="red", label="drop port")
axislegend(ax, position=:rt)

save("Ring_transmission.png", fig, px_per_unit = 2)