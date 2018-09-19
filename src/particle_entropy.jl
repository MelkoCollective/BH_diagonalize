"""
Calculate the particle entanglement entropy for a subset A, using the SVD.
"""
function particle_entropy(basis::AbstractSzbasis, Asize::Int, d::Vector{T}) where {T<:Number}
    basisA, basisB = particle_entropy_bases(basis, Asize)
    # Matrix to SVD
    Amatrix = zeros(T, length(basisA), length(basisB))

    for (i, braA) in enumerate(basisA)
        for (j, braB) in enumerate(basisB)
            bra = braA + braB
            bra in basis || continue

            norm = 1.0
            for k in 1:basis.K
                norm *= binomial(bra[k], braA[k])
            end

            Amatrix[i, j] = sqrt(norm) * d[serial_num(basis, bra)]
        end
    end

    S = svdvals(Amatrix) / sqrt(binomial(basis.N, basisA.N))
    err = abs(sum(S.^2) - 1.0)

    if err > 1e-12
        @warn("RDM eigenvalue error: $(err)")
    end

    -log(sum(S.^4))
end


"""
Generate the bases used for the particle entanglement entropy.
"""
function particle_entropy_bases(basis::Szbasis, Asize::Int)
    Bsize = basis.N - Asize
    basisA = Szbasis(basis.K, Asize)
    basisB = Szbasis(basis.K, Bsize)

    basisA, basisB
end

function particle_entropy_bases(basis::RestrictedSzbasis, Asize::Int)
    Bsize = basis.N - Asize
    if Asize <= basis.M
        basisA = Szbasis(basis.K, Asize)
    else
        basisA = RestrictedSzbasis(basis.K, Asize, basis.M)
    end
    if Bsize <= basis.M
        basisB = Szbasis(basis.K, Bsize)
    else
        basisB = RestrictedSzbasis(basis.K, Bsize, basis.M)
    end

    basisA, basisB
end
