module Pfapack77

using Libdl
using LinearAlgebra

global libp77 = C_NULL
global ilaenv = C_NULL

__init__() = begin
    global libp77 = dlopen("/opt/Pfapack77/lib/libpfapack.so")
    global ilaenv = dlsym(libp77, :ilaenv_)
end

include("ilaenv.jl")
include("sktdsm.jl")
include("sktrf.jl")

end
