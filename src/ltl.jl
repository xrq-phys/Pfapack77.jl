"Simplified caller. Allocated additional memory."
ltl_decomp!(A::Matrix{T}) where {T} = begin
    m, n = size(A)
    nb = nb_sytrd(T, "L", m)
    Work = zeros(T, nb * m)

    # Call.
    _, iPrm = sktrf!("L", "N", A, Work)

    Lcore = LowerTriangular(A[2:end, 1:end-1]);
    for i = 1:m-1
        Lcore[i, i] = 1.0;
    end
    L = [ 1. zeros(1, m-1);
          zeros(m-1) Lcore; ];
    Tv = [A[i+1, i] for i=1:m-1];

    (L, Tv, iPrm)
end

"Computes invert from LTL' decomposition."
ltl_inv(L, T, iPrm) = begin
    Li = inv(LowerTriangular(L))
    rP = invperm(iPrm)
    TLi = sktdsmx(T, Li)
    (Li' * TLi)[rP, rP]
end

"Computes invert by using LTL' decomposition."
ltl_inv!(A::Matrix{T}) where {T} = begin
    m, n = size(A)
    nb = nb_sytrd(T, "L", m)
    Work = zeros(T, nb * m)

    # Decompose.
    _, iPrm = sktrf!("L", "N", A, Work)

    # Invert right-lower block of L.
    Lcore = view(A, 2:m, 1:m-1)
    LAPACK.trtri!('L', 'U', Lcore)

    # Copy & clear clutter parts. L spans new memory.
    L = [ 1. zeros(1, m-1);
          zeros(m-1) Lcore; ];
    for j = 1:m
        L[j, j] = 1.0
        L[1:j-1, j] .= 0.0
    end

    # TODO: Use T w/o copying?
    Tv = [A[i+1, i] for i=1:m-1];
    sktdsmx!(A, Tv, L)
    rP = invperm(iPrm)

    A .= A[:, rP] # TODO: Include in sktdsmx.
    BLAS.trmm!('L', 'L', 'T', 'U', 1.0, L, A)
    A .= A[rP, :]
end

