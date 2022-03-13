ssktrf!(UpLo, Mode, A::Matrix{Float32}, iPiv::Vector{Cint}, Work::Array{Float32}) = ccall(dlsym(libp77, :ssktrf_), Cvoid,
    (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Float32}, Ptr{Cint}, Ptr{Cint}, Ptr{Float32}, Ptr{Cint}, Ptr{Cint}),
    UpLo, Mode, [size(A)[1]], A, [strides(A)[2]], iPiv, Work, [size(Work)[1]], [0])

dsktrf!(UpLo, Mode, A::Matrix{Float64}, iPiv::Vector{Cint}, Work::Array{Float64}) = ccall(dlsym(libp77, :dsktrf_), Cvoid,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}),
    UpLo, Mode, [size(A)[1]], A, [strides(A)[2]], iPiv, Work, [size(Work)[1]], [0])

csktrf!(UpLo, Mode, A::Matrix{ComplexF32}, iPiv::Vector{Cint}, Work::Array{ComplexF32}) = ccall(dlsym(libp77, :csktrf_), Cvoid,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}),
    UpLo, Mode, [size(A)[1]], A, [strides(A)[2]], iPiv, Work, [size(Work)[1]], [0])

zsktrf!(UpLo, Mode, A::Matrix{ComplexF64}, iPiv::Vector{Cint}, Work::Array{ComplexF64}) = ccall(dlsym(libp77, :zsktrf_), Cvoid,
        (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}, Ptr{Cint}),
    UpLo, Mode, [size(A)[1]], A, [strides(A)[2]], iPiv, Work, [size(Work)[1]], [0])

sktrf!(UpLo, Mode, A::Matrix{Float32}, iPiv::Vector{Cint}, Work::Array{Float32}) = ssktrf!(UpLo, Mode, A, iPiv, Work)
sktrf!(UpLo, Mode, A::Matrix{Float64}, iPiv::Vector{Cint}, Work::Array{Float64}) = dsktrf!(UpLo, Mode, A, iPiv, Work)
sktrf!(UpLo, Mode, A::Matrix{ComplexF32}, iPiv::Vector{Cint}, Work::Array{ComplexF32}) = csktrf!(UpLo, Mode, A, iPiv, Work)
sktrf!(UpLo, Mode, A::Matrix{ComplexF64}, iPiv::Vector{Cint}, Work::Array{ComplexF64}) = zsktrf!(UpLo, Mode, A, iPiv, Work)

sktrf!(UpLo, Mode, A::Matrix{T}, Work::Vector{T}) where {T} = begin
    m, n = size(A)
    m == n || throw(DimensionMismatch("Input matrix should be square."))
    iPiv = zeros(Cint, m)
    iPrm = [1:m...]
    sktrf!(UpLo, Mode, A, iPiv, Work)

    iRange = if lowercase(UpLo)[1] == 'l'
        1:1:m
    else
        m:-1:1
    end

    for i = iRange
        j = iPiv[i]
        index_0 = iPrm[i];
        iPrm[i] = iPrm[j];
        iPrm[j] = index_0;
    end
    (A, iPrm)
end

