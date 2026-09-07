#Testing macro
macro (add_test_executable_with_lib _name)
	add_executable(test_${_name} "${_name}.cpp")
	target_link_libraries(test_${_name} ${EMP-OT_LIBRARIES})
endmacro()

# ---------------------------------------------------------------------------
# Self-contained CTest launchers. Each one starts every party itself on a
# dedicated port block and fails if ANY party exits non-zero. They need no
# external `./run` wrapper (which is not in the repo), so `ctest` runs them
# directly. The field is argv[3]: fp61 | fp107 | f2k | z2k.
# ---------------------------------------------------------------------------

# 2-party drivers (vole_fp, cvole_fp): party 1 = ALICE, party 2 = BOB.
macro (add_2party_field_test _name _field _port)
	add_test(NAME ${_name}_${_field}
		COMMAND sh -c "B='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${_name}'; \"$B\" 1 ${_port} ${_field} & p1=$!; \"$B\" 2 ${_port} ${_field}; r2=$?; wait $p1; r1=$?; exit $((r1 | r2))"
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

# ncvole_fp: party 0 = server/committer, parties 1..2 = clients (n_client = 2).
macro (add_ncvole_field_test _field _port)
	add_test(NAME ncvole_fp_${_field}
		COMMAND sh -c "B='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_ncvole_fp'; \"$B\" 0 ${_port} ${_field} & s=$!; \"$B\" 1 ${_port} ${_field} & a=$!; \"$B\" 2 ${_port} ${_field} & b=$!; wait $a; r1=$?; wait $b; r2=$?; wait $s; r0=$?; exit $((r0 | r1 | r2))"
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

# mcvole_fp: n_party+1 processes, party ids 0..3 (king = party 3 = n_party).
macro (add_mcvole_field_test _field _port)
	add_test(NAME mcvole_fp_${_field}
		COMMAND sh -c "B='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_mcvole_fp'; \"$B\" 0 ${_port} ${_field} & a=$!; \"$B\" 1 ${_port} ${_field} & b=$!; \"$B\" 2 ${_port} ${_field} & c=$!; \"$B\" 3 ${_port} ${_field} & d=$!; wait $a; r0=$?; wait $b; r1=$?; wait $c; r2=$?; wait $d; r3=$?; exit $((r0 | r1 | r2 | r3))"
		WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

# Build the test executables (field chosen at run time via argv[3]).
add_test_executable_with_lib(vole_fp)
add_test_executable_with_lib(cvole_fp)
add_test_executable_with_lib(ncvole_fp)
add_test_executable_with_lib(mcvole_fp)
add_test_executable_with_lib(com_binding)
add_test_executable_with_lib(mpfss_chi)
add_test_executable_with_lib(vole_f2k)
add_test_executable_with_lib(vole_f2)

# Local (no network) MPFSS check-coefficient test: per-tree chi, agreement.
add_test(NAME mpfss_chi
	COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_mpfss_chi"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

# Local (no network) commitment checks: hiding sanity, shift-by-one, rank.
add_test(NAME com_binding
	COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_com_binding"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

# Register every driver over the default field (fp61) and the Z_{2^k} ring (z2k,
# which exercises the MozZ2karella ring consistency checks). Distinct port blocks
# (>=64 apart; mcvole's mesh uses ~24 ports) so the cases can run under `ctest -j`.
add_2party_field_test(vole_fp  fp61 23110)
add_2party_field_test(vole_fp  z2k  23210)
add_test(NAME vole_f2
	COMMAND sh -c "B='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vole_f2'; \"$B\" 1 24010 & p1=$!; \"$B\" 2 24010; r2=$?; wait $p1; r1=$?; exit $((r1 | r2))"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
add_test(NAME vole_f2k
	COMMAND sh -c "B='${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_vole_f2k'; \"$B\" 1 23910 & p1=$!; \"$B\" 2 23910; r2=$?; wait $p1; r1=$?; exit $((r1 | r2))"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
add_2party_field_test(cvole_fp fp61 23310)
add_2party_field_test(cvole_fp z2k  23410)
add_ncvole_field_test(fp61 23510)
add_ncvole_field_test(z2k  23610)
add_mcvole_field_test(fp61 23710)
add_mcvole_field_test(z2k  23810)
