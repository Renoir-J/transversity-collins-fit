using Libdl, LinearAlgebra, JSON
length(ARGS) >= 1 || error("usage: check_bridge.jl /absolute/libcollins_operator.so [report.json]")
library = Libdl.dlopen(abspath(ARGS[1]))
function_pointer = Libdl.dlsym(library, :cf_apply_operator)
U = [0.8 0.1 0.0; 0.1 0.7 0.2]
initial = [0.3 0.9; -0.2 0.5; 0.8 -0.1]
result = zeros(2, 2)
status = Ref{Cint}(-1)
ccall(function_pointer, Cvoid, (Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
      2, 3, 2, U, initial, result, status)
error_value = maximum(abs.(result - U * initial))
passed = status[] == 0 && error_value < 1e-14
report = Dict("passed" => passed, "scope" => "synthetic_bridge_ABI", "max_abs_error" => error_value,
              "HOPPET_transversity_operator_builder" => "NOT_IMPLEMENTED", "AD_route" => "Julia operator application")
if length(ARGS) > 1
    open(ARGS[2], "w") do stream
        JSON.print(stream, report, 2)
    end
end
println(JSON.json(report))
Libdl.dlclose(library)
passed || error("Fortran ABI validation failed")

