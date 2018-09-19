module BoseHubbardDiagonalize

using LinearAlgebra: svdvals
using SparseArrays: sparse

using JeszenszkiBasis

export
    BdryCond,
    OBC,
    PBC,

    sparse_hamiltonian,
    particle_entropy,
    spatial_entropy

"""
Boundary conditions.
"""
@enum BdryCond begin
    PBC
    OBC
end
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

include("sparse_hamiltonian.jl")
include("particle_entropy.jl")
include("spatial_entropy.jl")

end
