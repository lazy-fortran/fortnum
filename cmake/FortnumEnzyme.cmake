# FortnumEnzyme.cmake -- discover the Flang + LLVM + Enzyme toolchain.
#
# Sets FORTNUM_ENZYME_AVAILABLE (cache BOOL) after probing for the three
# binaries the Enzyme pipeline needs and the Enzyme LLVM pass plugin. When
# FORTNUM_ENABLE_ENZYME is ON and the plugin is missing, Enzyme tests skip
# unless FORTNUM_ENZYME_REQUIRED is ON, in which case configure hard-errors.
#
# This module owns discovery only. FortnumAddEnzymeTest.cmake consumes the
# result variables to build and register the actual derivative-check tests.
#
# Result/cache variables set here:
#   FORTNUM_FLANG_EXECUTABLE       Flang driver that emits LLVM IR and links.
#   FORTNUM_LLVM_OPT_EXECUTABLE    opt driver that runs the Enzyme pass.
#   FORTNUM_ENZYME_PLUGIN          Path to the Enzyme LLVM pass plugin (.so).
#   FORTNUM_LLVM_VERSION           Version string reported by opt.
#   FORTNUM_ENZYME_AVAILABLE       TRUE only when all of the above resolved.

include_guard(GLOBAL)

if(NOT FORTNUM_ENABLE_ENZYME)
    set(FORTNUM_ENZYME_AVAILABLE FALSE CACHE BOOL
        "Enzyme autodiff pipeline is usable" FORCE)
    return()
endif()

# Optional toolchain prefix. With no prefix, discovery uses PATH and the
# platform's standard library locations.
set(FORTNUM_LLVM_ROOT "$ENV{FORTNUM_LLVM_ROOT}" CACHE PATH
    "Optional LLVM toolchain prefix used by the Enzyme pipeline")

if(FORTNUM_LLVM_ROOT
        AND NOT DEFINED FORTNUM_FLANG_EXECUTABLE
        AND EXISTS "${FORTNUM_LLVM_ROOT}/bin/flang")
    set(FORTNUM_FLANG_EXECUTABLE "${FORTNUM_LLVM_ROOT}/bin/flang"
        CACHE FILEPATH "Flang driver used for LLVM IR emission and final link")
endif()
if(FORTNUM_LLVM_ROOT
        AND NOT DEFINED FORTNUM_CLANG_EXECUTABLE
        AND EXISTS "${FORTNUM_LLVM_ROOT}/bin/clang")
    set(FORTNUM_CLANG_EXECUTABLE "${FORTNUM_LLVM_ROOT}/bin/clang"
        CACHE FILEPATH "Clang driver used for Enzyme custom-rule helper IR")
endif()
if(FORTNUM_LLVM_ROOT
        AND NOT DEFINED FORTNUM_LLVM_OPT_EXECUTABLE
        AND EXISTS "${FORTNUM_LLVM_ROOT}/bin/opt")
    set(FORTNUM_LLVM_OPT_EXECUTABLE "${FORTNUM_LLVM_ROOT}/bin/opt"
        CACHE FILEPATH "LLVM opt driver used to run the Enzyme pass")
endif()

# --- Flang driver. Allow override via cache var. ---
find_program(FORTNUM_FLANG_EXECUTABLE
    NAMES flang-new flang
    DOC "Flang driver used for LLVM IR emission and final link")

# C helpers are used for Enzyme's custom derivative registration metadata.
find_program(FORTNUM_CLANG_EXECUTABLE
    NAMES clang clang-22 clang-21 clang-20 clang-19
    DOC "Clang driver used for Enzyme custom-rule helper IR")

# --- opt driver that loads the Enzyme pass. ---
find_program(FORTNUM_LLVM_OPT_EXECUTABLE
    NAMES opt opt-22 opt-21 opt-20 opt-19
    DOC "LLVM opt driver used to run the Enzyme pass")

# --- LLVM version, read from opt itself so it matches the pass plugin ABI. ---
set(FORTNUM_LLVM_VERSION "unknown")
if(FORTNUM_LLVM_OPT_EXECUTABLE)
    execute_process(
        COMMAND ${FORTNUM_LLVM_OPT_EXECUTABLE} --version
        OUTPUT_VARIABLE _fortnum_opt_version_out
        ERROR_VARIABLE _fortnum_opt_version_out
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "LLVM version ([0-9]+\\.[0-9]+\\.[0-9]+)"
        _fortnum_llvm_ver_match "${_fortnum_opt_version_out}")
    if(CMAKE_MATCH_1)
        set(FORTNUM_LLVM_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()

string(REGEX MATCH "^([0-9]+)" _fortnum_llvm_major "${FORTNUM_LLVM_VERSION}")
set(_fortnum_llvm_major "${CMAKE_MATCH_1}")

# --- Enzyme LLVM pass plugin. Search standard install prefixes; honor an
#     explicit FORTNUM_ENZYME_PLUGIN cache override unconditionally. ---
set(_fortnum_enzyme_plugin_names
    LLVMEnzyme-${_fortnum_llvm_major}.so
    LLVMEnzyme.so
    LLVMEnzyme-${_fortnum_llvm_major}.dylib
    LLVMEnzyme.dylib)

set(_fortnum_enzyme_search_paths
    /usr/lib
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64
    /opt/enzyme/lib)
if(FORTNUM_LLVM_ROOT)
    list(PREPEND _fortnum_enzyme_search_paths "${FORTNUM_LLVM_ROOT}/lib")
endif()

if(FORTNUM_LLVM_ROOT
        AND NOT DEFINED FORTNUM_ENZYME_PLUGIN
        AND EXISTS
            "${FORTNUM_LLVM_ROOT}/lib/LLVMEnzyme-${_fortnum_llvm_major}.so")
    set(FORTNUM_ENZYME_PLUGIN
        "${FORTNUM_LLVM_ROOT}/lib/LLVMEnzyme-${_fortnum_llvm_major}.so"
        CACHE FILEPATH
        "Enzyme LLVM pass plugin shared object (set manually to activate)")
endif()

find_library(FORTNUM_ENZYME_PLUGIN
    NAMES ${_fortnum_enzyme_plugin_names}
          LLVMEnzyme-${_fortnum_llvm_major}
          LLVMEnzyme
    PATHS ${_fortnum_enzyme_search_paths}
        ENV ENZYME_PLUGIN_DIR
    PATH_SUFFIXES enzyme Enzyme llvm/lib llvm-${_fortnum_llvm_major}/lib
    NO_DEFAULT_PATH
    DOC "Enzyme LLVM pass plugin shared object (set manually to activate)")

# --- Decide availability. ---
if(FORTNUM_FLANG_EXECUTABLE
        AND FORTNUM_LLVM_OPT_EXECUTABLE
        AND FORTNUM_ENZYME_PLUGIN
        AND EXISTS "${FORTNUM_ENZYME_PLUGIN}")
    set(FORTNUM_ENZYME_AVAILABLE TRUE CACHE BOOL
        "Enzyme autodiff pipeline is usable" FORCE)
else()
    set(FORTNUM_ENZYME_AVAILABLE FALSE CACHE BOOL
        "Enzyme autodiff pipeline is usable" FORCE)
endif()

# --- Report and gate. ---
if(FORTNUM_ENZYME_AVAILABLE)
    message(STATUS "fortnum: Enzyme available "
        "(flang=${FORTNUM_FLANG_EXECUTABLE}, opt=${FORTNUM_LLVM_OPT_EXECUTABLE}, "
        "llvm=${FORTNUM_LLVM_VERSION}, plugin=${FORTNUM_ENZYME_PLUGIN})")
else()
    set(_fortnum_enzyme_missing "")
    if(NOT FORTNUM_FLANG_EXECUTABLE)
        string(APPEND _fortnum_enzyme_missing " flang-new")
    endif()
    if(NOT FORTNUM_LLVM_OPT_EXECUTABLE)
        string(APPEND _fortnum_enzyme_missing " opt")
    endif()
    if(NOT FORTNUM_ENZYME_PLUGIN)
        string(APPEND _fortnum_enzyme_missing " Enzyme-plugin")
    endif()
    if(FORTNUM_ENZYME_REQUIRED)
        message(FATAL_ERROR
            "fortnum: FORTNUM_ENZYME_REQUIRED=ON but the Enzyme toolchain is "
            "incomplete; missing:${_fortnum_enzyme_missing}. Install the "
            "Enzyme pass plugin or point -DFORTNUM_ENZYME_PLUGIN=<path> at it. "
            "See docs/design/enzyme_toolchain.md.")
    else()
        message(STATUS
            "fortnum: Enzyme unavailable; missing:${_fortnum_enzyme_missing}. "
            "Enzyme derivative tests will be registered as skipped. "
            "See docs/design/enzyme_toolchain.md.")
    endif()
endif()
