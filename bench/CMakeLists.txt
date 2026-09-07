# Runtime benchmarks (separate from test/, which is for correctness only).
# Each is a 2-party program: bench_<x> <party> <port> [field] [threads] [size].
# Run both parties manually, e.g.:
#   ./bin/bench_cvole 1 12345 z2k 4 22 &
#   ./bin/bench_cvole 2 12345 z2k 4 22
macro (add_bench _name)
	add_executable(${_name} "${_name}.cpp")
	target_link_libraries(${_name} ${EMP-OT_LIBRARIES})
endmacro()

add_bench(bench_cvole)  # committed VOLE: throughput + A/B/shared split
add_bench(bench_vole)   # primal VOLE: throughput
add_bench(bench_mvole)  # primal n-party VOLE: n_party+1 processes, see the file header
