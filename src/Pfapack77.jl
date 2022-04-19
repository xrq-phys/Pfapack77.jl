module Pfapack77

using Libdl
using LinearAlgebra
using OpenBLAS32_jll
using Pfapack_jll: libpfapack

global libp77 = C_NULL
global ilaenv = C_NULL

__init__() = begin
    # Ensure 32-bit BLAS APIS are loaded.
    if (o -> o.interface != :lp64).(BLAS.lbt_get_config().loaded_libs) |> all
        BLAS.lbt_forward(OpenBLAS32_jll.libopenblas)
    end

    global libp77 = dlopen(libpfapack)
    global ilaenv = dlsym(libp77, :ilaenv_)
end

include("ilaenv.jl")
include("sktdsm.jl")
include("sktrf.jl")
include("ltl.jl")
include("utu.jl")

end
