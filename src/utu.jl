"Simplified caller. Allocated additional memory."
utu_decomp!(A::Matrix{T}) where {T} = begin
    m, n = size(A)
    nb = nb_sytrd(T, "U", m)
    Work = zeros(T, nb * m)

    # Call.
    _, iPrm = sktrf!("U", "N", A, Work)

    Ucore = UpperTriangular(A[1:end-1, 2:end]);
    for i = 1:m-1
        Ucore[i, i] = 1.0;
    end
    U = [ Ucore zeros(m-1);
          zeros(1, m-1) 1.; ];
    Tv = [A[i, i+1] for i=1:m-1];

    (U, Tv, iPrm)
end

"Computes invert from UTU' decomposition."
utu_inv(U, T, iPrm) = begin
    Ui = inv(UpperTriangular(U))
    rP = invperm(iPrm)
    TUi = sktdsmx(-T, Ui)
    (Ui' * TUi)[rP, rP]
end

"Computes invert by using UTU' decomposition."
utu_inv!(A::Matrix{T}) where {T} = begin
    m, n = size(A)
    nb = nb_sytrd(T, "U", m)
    Work = zeros(T, nb * m)

    # Decompose.
    _, iPrm = sktrf!("U", "N", A, Work)

    # Invert right-lower block of L.
    Ucore = view(A, 1:m-1, 2:m)
    LAPACK.trtri!('U', 'U', Ucore)

    # Copy & clear clutter parts. L spans new memory.
    U = [ Ucore zeros(m-1); 
          zeros(1, m-1) 1.; ];
    for i = 1:m
        U[i, i] = 1.0
        U[i, 1:i-1] .= 0.0
    end

    # TODO: Use T w/o copying?
    Tv = [A[i, i+1] for i=1:m-1];
    sktdsmx!(A, -Tv, U)
    rP = invperm(iPrm)

    A .= A[:, rP] # TODO: Include in sktdsmx.
    BLAS.trmm!('L', 'U', 'T', 'U', 1.0, U, A)
    A .= A[rP, :]
end

