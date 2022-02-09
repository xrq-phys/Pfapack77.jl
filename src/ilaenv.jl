nb_sytrd(::Type{Float32}, UpLo, n)    = ccall(ilaenv, Cint, (Ptr{Cint}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), [1], "SSYTRD", UpLo, [n], [-1], [-1], [-1])
nb_sytrd(::Type{Float64}, UpLo, n)    = ccall(ilaenv, Cint, (Ptr{Cint}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), [1], "DSYTRD", UpLo, [n], [-1], [-1], [-1])
nb_sytrd(::Type{ComplexF32}, UpLo, n) = ccall(ilaenv, Cint, (Ptr{Cint}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), [1], "CSYTRD", UpLo, [n], [-1], [-1], [-1])
nb_sytrd(::Type{ComplexF64}, UpLo, n) = ccall(ilaenv, Cint, (Ptr{Cint}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), [1], "ZSYTRD", UpLo, [n], [-1], [-1], [-1])
