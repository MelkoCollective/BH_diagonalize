# Operational entanglement peak for Bose-Hubbard chains in 1D.

using BoseHubbardDiagonalize

using ArgParse
using Arpack
using JeszenszkiBasis
using LsqFit

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "M"
        help = "number of sites"
        arg_type = Int
        required = true
    "N"
        help = "number of particles"
        arg_type = Int
        required = true
    "MA"
        help = "number of sites in region A"
        arg_type = Int
        required = true
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
end
add_arg_group(s, "boundary conditions")
@add_arg_table s begin
    "--pbc"
        help = "periodic boundary conditions (default)"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = PBC
        default = PBC
    "--obc"
        help = "open boundary conditions"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = OBC
        default = PBC
end
add_arg_group(s, "BH parameters")
@add_arg_table s begin
    "--u-min"
        metavar = "U"
        help = "minimum U"
        arg_type = Float64
        required = true
    "--u-max"
        metavar = "U"
        help = "maximum U"
        arg_type = Float64
        required = true
    "--ut-tol"
        metavar = "Ut"
        help = "required tolerance for U/t"
        arg_type = Float64
        default = 1e-3
    "--t"
        metavar = "t"
        help = "t value"
        arg_type = Float64
        default = 1.0
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Size of region A
const Asize = c[:MA]
# Output file
const output = c[:out]
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary]

if isnothing(site_max)
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

# Fitting model.
model(x, p) = p[1] .+ p[2]*x .+ p[3]*x.^2

open(output, "w") do f
    if isnothing(site_max)
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")
    end
    write(f, "# U/t E0/t Eop(l=$(Asize))\n")

    Us = Float64[]
    S2s = Float64[]
    U = c[:u_min]

    while true
        H = sparse_hamiltonian(basis, c[:t], U, boundary=boundary)
        d = eigs(H, nev=1, which=:SR)
        wf = vec(d[2])
        _, s2_operational = spatial_entropy(basis, Asize, wf)

        write(f, "$(U/c[:t]) $(d[1][1]/c[:t]) $(s2_operational)\n")
        flush(f)

        push!(Us, U)
        push!(S2s, s2_operational)

        if length(Us) == 1
            U = c[:u_max]
        elseif length(Us) == 2
            U = 0.5 * (Us[1] + Us[2])
        else
            # Find the maximum and the two closest points around it.
            U_max_idx = findmax(S2s)[2]
            U_max = Us[U_max_idx]
            if isempty(Us[Us .< U_max]) || isempty(Us[Us .> U_max])
                error("No window")
            end

            U_left, U_left_idx = findmax([U < U_max ? U : -Inf for U in Us])
            U_right, U_right_idx = findmin([U > U_max ? U : Inf for U in Us])

            # If we have tight enough bounds, we're done.
            bound = max(U_max - U_left, U_right - U_max) / c[:t]
            println("Within $(bound) at $(U/c[:t]), $(s2_operational) (max at $(U_max/c[:t]), $(S2s[U_max_idx]))")
            if bound <= c[:ut_tol]
                break
            end

            # We need to get closer, so try a parabolic fit based on the 3
            # points we've selected.
            xs = [U_left, U_max, U_right]
            ys = S2s[[U_left_idx, U_max_idx, U_right_idx]]
            fit = curve_fit(model, xs, ys, [maximum(S2s), 1.0, -1.0])
            U = -0.5 * fit.param[2] / fit.param[3]

            # If the fit isn't reasonable, we don't use it at all.
            if !(U_left < U < U_right)
                @info("Bad fit: $(U/c[:t]), bisecting")
                if S2s[U_left_idx] > S2s[U_right_idx]
                    U = 0.5 * (U_left + U_max)
                else
                    U = 0.5 * (U_max + U_right)
                end
            end
        end
    end
end
