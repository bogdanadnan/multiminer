cmake_minimum_required (VERSION 2.8)
project (a.multiminer C ASM CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

# Windows build additional configuration
# set(DEPRECATION_FLAGS "-Xlinker --allow-multiple-definition -Wno-error=deprecated")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEPRECATION_FLAGS}")
# set(C_DEP_FLAGS "-Wno-pointer-sign -Wno-pointer-to-int-cast")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEP_FLAGS}")
# set(NO_CUDA TRUE)

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "./")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(FindJansson.cmake)

set(SOURCE
        cpu-miner.c
        util.c
        uint256.cpp
        api.c
        sysinfos.c
        algo-gate-api.c
        compat/winansi.c
        crypto/oaes_lib.c
        crypto/c_keccak.c
        crypto/c_groestl.c
        crypto/c_blake256.c
        crypto/c_jh.c
        crypto/c_skein.c
        crypto/hash.c
        crypto/aesb.c
        crypto/magimath.cpp
        algo/argon2/argon2a/argon2a.c
        algo/argon2/argon2a/ar2/argon2.c
        algo/argon2/argon2a/ar2/opt.c
        algo/argon2/argon2a/ar2/cores.c
        algo/argon2/argon2a/ar2/ar2-scrypt-jane.c
        algo/argon2/argon2a/ar2/blake2b.c
        algo/argon2/argon2d/argon2d-gate.c
        algo/argon2/argon2d/blake2/blake2b.c
        algo/argon2/argon2d/argon2d/argon2.c
        algo/argon2/argon2d/argon2d/core.c
        algo/argon2/argon2d/argon2d/opt.c
        algo/argon2/argon2d/argon2d/thread.c
        algo/argon2/argon2d/argon2d/encoding.c
        algo/blake/sph_blake.c
        algo/blake/blake-hash-4way.c
        algo/blake/blake-gate.c
        algo/blake/blake.c
        algo/blake/blake-4way.c
        algo/blake/sph_blake2b.c
        algo/blake/blake2b.c
        algo/blake/sph-blake2s.c
        algo/blake/blake2s-hash-4way.c
        algo/blake/blake2s.c
        algo/blake/blake2s-gate.c
        algo/blake/blake2s-4way.c
        algo/blake/blakecoin-gate.c
        algo/blake/mod_blakecoin.c
        algo/blake/blakecoin.c
        algo/blake/blakecoin-4way.c
        algo/blake/decred-gate.c
        algo/blake/decred.c
        algo/blake/decred-4way.c
        algo/blake/pentablake-gate.c
        algo/blake/pentablake-4way.c
        algo/blake/pentablake.c
        algo/bmw/sph_bmw.c
        algo/bmw/bmw-hash-4way.c
        algo/bmw/bmw256.c
        algo/cryptonight/cryptolight.c
#        algo/cryptonight/cryptonight-common.c
#        algo/cryptonight/cryptonight-aesni.c
#        algo/cryptonight/cryptonight.c
        algo/cubehash/sph_cubehash.c
        algo/cubehash/sse2/cubehash_sse2.c
        algo/cubehash/cube-hash-2way.c
        algo/echo/sph_echo.c
        algo/echo/aes_ni/hash.c
        algo/gost/sph_gost.c
        algo/groestl/sph_groestl.c
        algo/groestl/groestl.c
        algo/groestl/myrgr-gate.c
        algo/groestl/myrgr-4way.c
        algo/groestl/myr-groestl.c
        algo/groestl/aes_ni/hash-groestl.c
        algo/groestl/aes_ni/hash-groestl256.c
        algo/fugue/sph_fugue.c
        algo/hamsi/sph_hamsi.c
        algo/hamsi/hamsi-hash-4way.c
        algo/haval/haval.c
        algo/haval/haval-hash-4way.c
        algo/heavy/sph_hefty1.c
        algo/heavy/heavy.c
        algo/heavy/bastion.c
        algo/hodl/aes.c
        algo/hodl/hodl-gate.c
        algo/hodl/hodl-wolf.c
        algo/hodl/sha512_avx.c
        algo/hodl/sha512_avx2.c
        algo/jh/sph_jh.c
        algo/jh/jh-hash-4way.c
        algo/jh/jha-gate.c
        algo/jh/jha-4way.c
        algo/jh/jha.c
        algo/keccak/sph_keccak.c
        algo/keccak/keccak.c
        algo/keccak/keccak-hash-4way.c
        algo/keccak/keccak-4way.c
        algo/keccak/keccak-gate.c
        algo/keccak/sse2/keccak.c
        algo/luffa/sph_luffa.c
        algo/luffa/luffa.c
        algo/luffa/luffa_for_sse2.c
        algo/luffa/luffa-hash-2way.c
        algo/lyra2/lyra2.c
        algo/lyra2/sponge.c
        algo/lyra2/lyra2rev2-gate.c
        algo/lyra2/lyra2rev2.c
        algo/lyra2/lyra2rev2-4way.c
        algo/lyra2/lyra2re.c
        algo/lyra2/lyra2z-gate.c
        algo/lyra2/lyra2z.c
        algo/lyra2/lyra2z-4way.c
        algo/lyra2/lyra2z330.c
        algo/lyra2/lyra2h-gate.c
        algo/lyra2/lyra2h.c
        algo/lyra2/lyra2h-4way.c
        algo/lyra2/allium-gate.c
        algo/lyra2/allium-4way.c
        algo/lyra2/allium.c
        algo/m7m.c
        algo/neoscrypt/neoscrypt.c
        algo/nist5/nist5-gate.c
        algo/nist5/nist5-4way.c
        algo/nist5/nist5.c
        algo/nist5/zr5.c
        algo/pluck.c
        algo/quark/quark-gate.c
        algo/quark/quark.c
        algo/quark/quark-4way.c
        algo/quark/anime-gate.c
        algo/quark/anime.c
        algo/quark/anime-4way.c
        algo/qubit/qubit-gate.c
        algo/qubit/qubit.c
        algo/qubit/qubit-2way.c
        algo/qubit/deep-gate.c
        algo/qubit/deep-2way.c
        algo/qubit/deep.c
        algo/ripemd/sph_ripemd.c
        algo/ripemd/ripemd-hash-4way.c
        algo/ripemd/lbry-gate.c
        algo/ripemd/lbry.c
        algo/ripemd/lbry-4way.c
        algo/scrypt.c
        algo/scryptjane/scrypt-jane.c
        algo/sha/sph_sha2.c
        algo/sha/sph_sha2big.c
        algo/sha/sha2-hash-4way.c
        algo/sha/sha2.c
        algo/sha/sha256t-gate.c
        algo/sha/sha256t-4way.c
        algo/sha/sha256t.c
        algo/shabal/sph_shabal.c
        algo/shabal/shabal-hash-4way.c
        algo/shavite/sph_shavite.c
        algo/shavite/sph-shavite-aesni.c
        algo/shavite/shavite.c
        algo/simd/sph_simd.c
        algo/simd/nist.c
        algo/simd/vector.c
        algo/simd/simd-hash-2way.c
        algo/skein/sph_skein.c
        algo/skein/skein-hash-4way.c
        algo/skein/skein.c
        algo/skein/skein-4way.c
        algo/skein/skein-gate.c
        algo/skein/skein2.c
        algo/skein/skein2-4way.c
        algo/skein/skein2-gate.c
        algo/sm3/sm3.c
        algo/sm3/sm3-hash-4way.c
        algo/tiger/sph_tiger.c
        algo/whirlpool/sph_whirlpool.c
        algo/whirlpool/whirlpool-hash-4way.c
        algo/whirlpool/whirlpool-gate.c
        algo/whirlpool/whirlpool-4way.c
        algo/whirlpool/whirlpool.c
        algo/whirlpool/whirlpoolx.c
        algo/x11/x11-gate.c
        algo/x11/x11.c
        algo/x11/x11-4way.c
        algo/x11/x11gost-gate.c
        algo/x11/x11gost.c
        algo/x11/x11gost-4way.c
        algo/x11/c11-gate.c
        algo/x11/c11.c
        algo/x11/c11-4way.c
        algo/x11/tribus-gate.c
        algo/x11/tribus.c
        algo/x11/tribus-4way.c
        algo/x11/timetravel-gate.c
        algo/x11/timetravel.c
        algo/x11/timetravel-4way.c
        algo/x11/timetravel10-gate.c
        algo/x11/timetravel10.c
        algo/x11/timetravel10-4way.c
        algo/x11/fresh.c
        algo/x11/x11evo.c
        algo/x11/x11evo-4way.c
        algo/x11/x11evo-gate.c
        algo/x12/x12-gate.c
        algo/x12/x12.c
        algo/x12/x12-4way.c
        algo/x13/x13-gate.c
        algo/x13/x13.c
        algo/x13/x13-4way.c
        algo/x13/x13sm3-gate.c
        algo/x13/x13sm3.c
        algo/x13/x13sm3-4way.c
        algo/x13/phi1612-gate.c
        algo/x13/phi1612.c
        algo/x13/phi1612-4way.c
        algo/x13/skunk-gate.c
        algo/x13/skunk-4way.c
        algo/x13/skunk.c
        algo/x13/drop.c
        algo/x14/x14-gate.c
        algo/x14/x14.c
        algo/x14/x14-4way.c
        algo/x14/veltor-gate.c
        algo/x14/veltor.c
        algo/x14/veltor-4way.c
        algo/x14/polytimos-gate.c
        algo/x14/polytimos.c
        algo/x14/polytimos-4way.c
#        algo/x14/axiom.c
        algo/x15/x15-gate.c
        algo/x15/x15.c
        algo/x15/x15-4way.c
        algo/x17/x17-gate.c
        algo/x17/x17.c
        algo/x17/x17-4way.c
        algo/x17/xevan-gate.c
        algo/x17/xevan.c
        algo/x17/xevan-4way.c
        algo/x17/x16r-gate.c
        algo/x17/x16r.c
        algo/x17/x16r-4way.c
        algo/x17/hmq1725.c
        algo/yescrypt/yescrypt.c
        algo/yescrypt/sha256_Y.c
        algo/yescrypt/yescrypt-best.c
        asm/sha2-x64.S asm/scrypt-x64.S asm/aesb-x64.S asm/neoscrypt_asm.S
        )

include(CheckIncludeFiles)
include(CheckFunctionExists)
include(CheckPrototypeDefinition)

check_include_files ( alloca.h HAVE_ALLOCA_H )
check_include_files ( inttypes.h HAVE_INTTYPES_H )
check_include_files ( memory.h HAVE_MEMORY_H )
check_include_files ( stdint.h HAVE_STDINT_H )
check_include_files ( stdlib.h HAVE_STDLIB_H )
check_include_files ( strings.h HAVE_STRINGS_H )
check_include_files ( string.h HAVE_STRING_H )
check_include_files ( syslog.h HAVE_SYSLOG_H )
check_include_files ( sys/endian.h HAVE_SYS_ENDIAN_H )
check_include_files ( sys/param.h HAVE_SYS_PARAM_H )
check_include_files ( sys/stat.h HAVE_SYS_STAT_H )
check_include_files ( sys/sysctl.h HAVE_SYS_SYSCTL_H )
check_include_files ( sys/types.h HAVE_SYS_TYPES_H )
check_include_files ( unistd.h HAVE_UNISTD_H )

check_function_exists ( alloca HAVE_ALLOCA )
check_function_exists ( getopt_long HAVE_GETOPT_LONG )

check_prototype_definition(be32dec "uint32_t be32dec(const void *)" "0" "sys/endian.h" HAVE_DECL_BE32DEC)
check_prototype_definition(be32enc "void be32enc(void *, uint32_t)" ""  "sys/endian.h" HAVE_DECL_BE32ENC)
check_prototype_definition(le32dec "uint32_t le32dec(const void *)" "0" "sys/endian.h" HAVE_DECL_LE32DEC)
check_prototype_definition(le32enc "void le32enc(void *, uint32_t)" "" "sys/endian.h" HAVE_DECL_LE32ENC)
check_prototype_definition(le16dec "uint16_t le16dec(const void *)" "0" "sys/endian.h" HAVE_DECL_LE16DEC)
check_prototype_definition(le16enc "void le16enc(void *, uint16_t)" "" "sys/endian.h" HAVE_DECL_LE16ENC)

find_package (Threads)
find_package(CURL)
find_package(Jansson)
find_package(OpenSSL)
find_package(GMP)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cpuminer-config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/cpuminer-config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(algo/argon2/argon2d)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CURL_INCLUDE_DIR})
include_directories(${JANSSON_INCLUDE_DIRS})
include_directories(${OPENSSL_INCLUDE_DIR})
include_directories(${GMP_INCLUDES})

add_executable(a.multiminer ${SOURCE})
target_link_libraries(a.multiminer mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
IF (WIN32)
	target_link_libraries(a.multiminer -lws2_32)
ENDIF()
# Windows build works with O1 instead of O3
target_compile_options(a.multiminer PRIVATE -O3 -mtune=native -march=native -I.)

IF (BUILD_ALL_ARCH)
	add_executable(a.multiminer-sse2 ${SOURCE})
	target_link_libraries(a.multiminer-sse2 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-sse2 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-sse2 PRIVATE -O3 -msse2 -I.)

	add_executable(a.multiminer-ssse3 ${SOURCE})
	target_link_libraries(a.multiminer-ssse3 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-ssse3 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-ssse3 PRIVATE -O3 -march=core2 -I.)

	add_executable(a.multiminer-sse42 ${SOURCE})
	target_link_libraries(a.multiminer-sse42 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-sse42 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-sse42 PRIVATE -O3 -march=corei7 -I.)

	add_executable(a.multiminer-aes-sse42 ${SOURCE})
	target_link_libraries(a.multiminer-aes-sse42 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-aes-sse42 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-aes-sse42 PRIVATE -O3 -maes -msse4.2 -I.)

	add_executable(a.multiminer-aes-avx ${SOURCE})
	target_link_libraries(a.multiminer-aes-avx mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-aes-avx -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-aes-avx PRIVATE -O3 -march=corei7-avx -I.)

	add_executable(a.multiminer-avx2 ${SOURCE})
	target_link_libraries(a.multiminer-avx2 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-avx2 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-avx2 PRIVATE -O3 -march=core-avx2 -I.)

	add_executable(a.multiminer-avx512 ${SOURCE})
	target_link_libraries(a.multiminer-avx512 mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-avx512 -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-avx512 PRIVATE -O3 -march=skylake-avx512 -fno-asynchronous-unwind-tables -I.)

	add_executable(a.multiminer-avx2-sha ${SOURCE})
	target_link_libraries(a.multiminer-avx2-sha mm_gpu_gate ${CMAKE_THREAD_LIBS_INIT} ${CURL_LIBRARIES} ${JANSSON_LIBRARIES} ${OPENSSL_CRYPTO_LIBRARY} ${OPENSSL_SSL_LIBRARY} ${GMP_LIBRARIES})
	IF (WIN32)
		target_link_libraries(a.multiminer-avx2-sha -lws2_32)
	ENDIF()
	target_compile_options(a.multiminer-avx2-sha PRIVATE -O3 -march=core-avx2 -I.)
ENDIF()

add_custom_command(
        TARGET a.multiminer POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/algo/argon2/argon2d/argon2-gpu/data/kernels/argon2_kernel.cl
                ${CMAKE_CURRENT_BINARY_DIR}/data/kernels/argon2_kernel.cl)

