module Pfapack77

using Libdl
using LinearAlgebra
if Int == Int64
    using OpenBLAS32
end
using Pfapack_jll: libpfapack

global libp77 = C_NULL
global ilaenv = C_NULL

__init__() = begin
    global libp77 = dlopen(libpfapack)
    global ilaenv = dlsym(libp77, :ilaenv_)
end

include("ilaenv.jl")
include("sktdsm.jl")
include("sktrf.jl")
include("ltl.jl")
include("utu.jl")

end
