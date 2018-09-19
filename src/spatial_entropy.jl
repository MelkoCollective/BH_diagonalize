"""
Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""
function spatial_entropy(basis::AbstractSzbasis, A, d::Vector{T}) where {T<:Number}
    B = setdiff(1:basis.K, A)

    isempty(A) && return 0.0, 0.0
    isempty(B) && return 0.0, 0.0

    # Matrices to SVD
    Amatrices = []
    for i in 0:basis.N
        DimA = num_vectors(basis, i, length(A))
        DimB = num_vectors(basis, basis.N-i, length(B))

        push!(Amatrices, zeros(T, DimA, DimB))
    end

    norms = zeros(Float64, basis.N+1)

    for (i, bra) in enumerate(basis)
        braA = view(bra, A)
        braB = view(bra, B)

        row = serial_num(basis, length(A), sum(braA), braA)
        col = serial_num(basis, length(B), sum(braB), braB)

        Amatrices[1 + sum(braA)][row, col] = d[i]
        norms[1 + sum(braA)] += abs2(d[i])
    end

    norm_err = abs(sum(norms) - 1.0)

    if norm_err > 1e-12
        @warn("norm error: $(norm_err)")
    end

    Ss_raw = [svdvals(Amatrix) for Amatrix in Amatrices]

    # Spatial.
    S_sp = vcat(Ss_raw...)
    err_sp = abs(sum(S_sp.^2) - 1.0)

    if err_sp > 1e-12
        @warn("RDM eigenvalue error: $(err_sp)")
    end

    S2_sp = -log(sum(S_sp.^4))

    # Operational.
    Ss_op = [S / sqrt(n) for (S, n) in zip(Ss_raw, norms)]
    errs_op = [abs(sum(S.^2) - 1.0) for S in Ss_op if !isempty(S)]

    if any(errs_op .> 1e-12)
        @warn("RDM eigenvalue error $(maximum(errs_op))")
    end

    S2_op = 0.0
    for (S, n) in zip(Ss_op, norms)
        isempty(S) && continue
        S2_op -= log(sum(S.^4)) * n
    end

    S2_sp, S2_op
end

spatial_entropy(basis::AbstractSzbasis, Asize::Int, d::Vector{T}) where {T<:Number} = spatial_entropy(basis, 1:Asize, d)
