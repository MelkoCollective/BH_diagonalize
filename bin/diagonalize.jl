# Renyi entanglement entropy of Bose-Hubbard chains in 1D.

using BoseHubbardDiagonalize

using ArgParse
using Arpack
using JeszenszkiBasis
using ProgressMeter: Progress, ProgressWrapper

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
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
end
add_arg_group(s, "output settings")
@add_arg_table s begin
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--no-progress"
        help = "hide progress bar"
        action = :store_true
    "--verbose", "-v"
        help = "show extra information"
        action = :store_true
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
        default = 1.0
    "--u-max"
        metavar = "U"
        help = "maximum U"
        arg_type = Float64
        default = 20.0
    "--u-step"
        metavar = "U"
        help = "U step"
        arg_type = Float64
    "--u-num"
        metavar = "N"
        help = "number of U"
        arg_type = Int
    "--u-log"
        help = "use logarithmic scale for U"
        action = :store_true
    "--t"
        metavar = "t"
        help = "t value"
        arg_type = Float64
        default = 1.0
end
add_arg_group(s, "entanglement entropy")
@add_arg_table s begin
    "--ee"
        metavar = "XA"
        help = "compute all EEs with partition size XA"
        arg_type = Int
        required = true
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Output file
const output = c[:out]
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary]
# Size of region A
const Asize = c[:ee]

const progress_output = c[:no_progress] ? devnull : stderr

try
    global U_range = range(c[:u_min]; stop=c[:u_max], length=c[:u_num], step=c[:u_step])
catch ex
    ex isa ArgumentError || rethrow()
    println("--u-step and --u-num may not both be supplied")
    exit(1)
end

c[:u_log] && (U_range = 10 .^ U_range)

if isnothing(site_max)
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

# Initial vector for diagonalization
const v0 = ones(Float64, length(basis))

# Diagonalization diagnostics
const niters = zeros(Int, length(U_range))
const nmults = zeros(Int, length(U_range))

open(output, "w") do f
    if isnothing(site_max)
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")
    end
    write(f, "# U/t E0/t S2(n=$(Asize)) S2(l=$(Asize)) Eop(l=$(Asize))\n")

    meter = Progress(length(U_range), output=progress_output)
    for (i, U) in ProgressWrapper(enumerate(U_range), meter)
        # Create the Hamiltonian
        H = sparse_hamiltonian(basis, c[:t], U, boundary=boundary)

        # Perform the Lanczos diagonalization to obtain the lowest eigenvector
        d = eigs(H, nev=1, which=:SR, v0=v0)
        E0 = d[1][1]
        wf = vec(d[2])
        d[3] == 1 || @warn("Diagonalization did not converge")
        niters[i] = d[4]
        nmults[i] = d[5]

        # Use the current ground state for the next diagonalization
        v0 .= wf

        # Calculate the second Renyi entropy
        s2_particle = particle_entropy(basis, Asize, wf)
        s2_spatial, s2_operational = spatial_entropy(basis, Asize, wf)

        write(f, "$(U/c[:t]) $(E0/c[:t]) $(s2_particle) $(s2_spatial) $(s2_operational)\n")
        flush(f)
    end
end

if c[:verbose]
    # Output diagnostics
    println("niter min/med/max = $(minimum(niters))/$(ceil(Int, median(niters)))/$(maximum(niters))")
    println("nmult min/med/max = $(minimum(nmults))/$(ceil(Int, median(nmults)))/$(maximum(nmults))")
end
