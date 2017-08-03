# Qubit-like entanglement in Bose-Hubbard chains in 1D.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "../src"))
using BoseHubbardDiagonalize

using ArgParse
using JeszenszkiBasis
using Qutilities

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
    "QA1"
        help = "first site of qubit A"
        arg_type = Int
        required = true
    "QA2"
        help = "second site of qubit A"
        arg_type = Int
        required = true
    "QB1"
        help = "first site of qubit B"
        arg_type = Int
        required = true
    "QB2"
        help = "second site of qubit B"
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

qubit_sites_A = [c[:QA1], c[:QA2]]
qubit_sites_B = [c[:QB1], c[:QB2]]
qubit_sites = [qubit_sites_A ; qubit_sites_B]
if any(qubit_sites .< 1) || any(qubit_sites .> M)
    println("Invalid qubit site")
    exit(1)
end

const bulk = setdiff(1:M, qubit_sites)

if c[:u_log] && c[:u_num] === nothing
    println("--u-log must be used with --u-num")
    exit(1)
end

if c[:u_step] === nothing
    if c[:u_num] === nothing
        U_range = c[:u_min]:0.5:c[:u_max]
    else
        if c[:u_log]
            U_range = logspace(c[:u_min], c[:u_max], c[:u_num])
        else
            U_range = linspace(c[:u_min], c[:u_max], c[:u_num])
        end
    end
else
    if c[:u_num] === nothing
        U_range = c[:u_min]:c[:u_step]:c[:u_max]
    else
        println("--u-step and --u-num may not both be supplied")
        exit(1)
    end
end

if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

println(join([i in qubit_sites_A ? "A" : (i in qubit_sites_B ? "B" : "C") for i in 1:M], " "))

open(output, "w") do f
    if site_max === nothing
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")
    end
    write(f, "# U/t E0/t P I2 Ivn EF EN pur_AB pur_A\n")

    for U in U_range
        H = sparse_hamiltonian(basis, c[:t], U, boundary=boundary)
        d = eigs(H, nev=1, which=:SR)
        wf = vec(d[2])

        # Weight of the projected subspace.
        P = 0.0
        # Reduced density matrix.
        rho_AB = zeros(Float64, 4, 4)
        for (i, v1) in enumerate(basis)
            sum(v1[qubit_sites_A]) == 1 && sum(v1[qubit_sites_B]) == 1 || continue
            row = 1 + 2 * Int(v1[qubit_sites_A] == [1, 0]) + Int(v1[qubit_sites_B] == [1, 0])

            P += wf[i]^2

            for (j, v2) in enumerate(basis)
                v1[bulk] == v2[bulk] || continue
                sum(v2[qubit_sites_A]) == 1 && sum(v2[qubit_sites_B]) == 1 || continue
                col = 1 + 2 * Int(v2[qubit_sites_A] == [1, 0]) + Int(v2[qubit_sites_B] == [1, 0])

                rho_AB[row, col] += wf[i] * wf[j]
            end
        end
        if P == 0.0
            error("Empty subspace")
        end
        rho_AB /= P

        err_trace = abs(trace(rho_AB) - 1.0)
        if err_trace > 1e-15
            warn("Bad trace: $(err_trace)")
        end

        # Very reduced density matrix.
        rho_A = ptrace(rho_AB)

        # Various entanglement measures.
        I2 = mutinf(rho_AB, S_renyi)
        Ivn = mutinf(rho_AB)
        EF = formation(rho_AB)
        EN = negativity(rho_AB)
        pur_AB = purity(rho_AB)
        pur_A = purity(rho_A)

        write(f, "$(U/c[:t]) $(d[1][1]/c[:t]) $(P) $(I2) $(Ivn) $(EF) $(EN) $(pur_AB) $(pur_A)\n")
        flush(f)
    end
end
