#[======================================================[.rst

Adds the following targets::

    pybind11::pybind11 - link to Python headers and pybind11::headers
    pybind11::module - Adds module links
    pybind11::embed - Adds embed links
    pybind11::lto - Link time optimizations (only if CMAKE_INTERPROCEDURAL_OPTIMIZATION is not set)
    pybind11::thin_lto - Link time optimizations (only if CMAKE_INTERPROCEDURAL_OPTIMIZATION is not set)
    pybind11::python_link_helper - Adds link to Python libraries
    pybind11::windows_extras - MSVC bigobj and mp for building multithreaded
    pybind11::opt_size - avoid optimizations that increase code size

Adds the following functions::

    pybind11_strip(target) - strip target after building on linux/macOS
    pybind11_find_import(module) - See if a module is installed.

#]======================================================]

# If we are in subdirectory mode, all IMPORTED targets must be GLOBAL. If we
# are in CONFIG mode, they should be "normal" targets instead.
# In CMake 3.11+ you can promote a target to global after you create it,
# which might be simpler than this check.
get_property(
  is_config
  TARGET pybind11::headers
  PROPERTY IMPORTED)

if(NOT is_config)
  include_guard(GLOBAL)
  set(optional_global GLOBAL)
else()
  include_guard(DIRECTORY)
  set(optional_global "")
endif()

# If not run in Python mode, we still would like this to at least
# include pybind11's include directory:
set(pybind11_INCLUDE_DIRS
    "${pybind11_INCLUDE_DIR}"
    CACHE INTERNAL "Include directory for pybind11 (Python not requested)")

# CMP0190 prohibits calling FindPython with both Interpreter and Development components
# when cross-compiling, unless the CMAKE_CROSSCOMPILING_EMULATOR variable is defined.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.1")
  cmake_policy(GET CMP0190 _pybind11_cmp0190)
  if(_pybind11_cmp0190 STREQUAL "NEW")
    set(PYBIND11_USE_CROSSCOMPILING "ON")
  endif()
endif()

if(CMAKE_CROSSCOMPILING
   AND PYBIND11_USE_CROSSCOMPILING
   AND NOT DEFINED CMAKE_CROSSCOMPILING_EMULATOR)
  set(_PYBIND11_CROSSCOMPILING
      ON
      CACHE INTERNAL "")
else()
  set(_PYBIND11_CROSSCOMPILING
      OFF
      CACHE INTERNAL "")
endif()

# --------------------- Shared targets ----------------------------

# Build an interface library target:
add_library(pybind11::pybind11 IMPORTED INTERFACE ${optional_global})
set_property(
  TARGET pybind11::pybind11
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES pybind11::headers)

# Build a module target:
add_library(pybind11::module IMPORTED INTERFACE ${optional_global})
set_property(
  TARGET pybind11::module
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES pybind11::pybind11)

# Build an embed library target:
add_library(pybind11::embed IMPORTED INTERFACE ${optional_global})
set_property(
  TARGET pybind11::embed
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES pybind11::pybind11)

# -------------- emscripten requires exceptions enabled -------------
# _pybind11_no_exceptions is a private mechanism to disable this addition.
# Please open an issue if you need to use it; it will be removed if no one
# needs it.
if(CMAKE_SYSTEM_NAME MATCHES Emscripten AND NOT _pybind11_no_exceptions)
  if(is_config)
    set(_tmp_config_target pybind11::pybind11_headers)
  else()
    set(_tmp_config_target pybind11_headers)
  endif()

  set_property(
    TARGET ${_tmp_config_target}
    APPEND
    PROPERTY INTERFACE_LINK_OPTIONS -fexceptions)
  set_property(
    TARGET ${_tmp_config_target}
    APPEND
    PROPERTY INTERFACE_COMPILE_OPTIONS -fexceptions)
  unset(_tmp_config_target)
endif()

# --------------------------- link helper ---------------------------

add_library(pybind11::python_link_helper IMPORTED INTERFACE ${optional_global})

set_property(
  TARGET pybind11::python_link_helper
  APPEND
  PROPERTY INTERFACE_LINK_OPTIONS "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>")

# ------------------------ Windows extras -------------------------

add_library(pybind11::windows_extras IMPORTED INTERFACE ${optional_global})

if(MSVC) # That's also clang-cl
  # /bigobj is needed for bigger binding projects due to the limit to 64k
  # addressable sections
  set_property(
    TARGET pybind11::windows_extras
    APPEND
    PROPERTY INTERFACE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)

  # /MP enables multithreaded builds (relevant when there are many files) for MSVC
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC") # no Clang no Intel
    # Only set these options for C++ files.  This is important so that, for
    # instance, projects that include other types of source files like CUDA
    # .cu files don't get these options propagated to nvcc since that would
    # cause the build to fail.
    set_property(
      TARGET pybind11::windows_extras
      APPEND
      PROPERTY INTERFACE_COMPILE_OPTIONS $<$<NOT:$<CONFIG:Debug>>:$<$<COMPILE_LANGUAGE:CXX>:/MP>>)
  endif()
endif()

# ----------------------- Optimize binary size --------------------------

add_library(pybind11::opt_size IMPORTED INTERFACE ${optional_global})

if(MSVC)
  set(PYBIND11_OPT_SIZE /Os)
else()
  set(PYBIND11_OPT_SIZE -Os)
endif()

set_property(
  TARGET pybind11::opt_size
  APPEND
  PROPERTY INTERFACE_COMPILE_OPTIONS $<$<CONFIG:Release>:${PYBIND11_OPT_SIZE}>
           $<$<CONFIG:MinSizeRel>:${PYBIND11_OPT_SIZE}>
           $<$<CONFIG:RelWithDebInfo>:${PYBIND11_OPT_SIZE}>)

# ----------------------- Legacy option --------------------------

# Warn or error if old variable name used
if(PYBIND11_CPP_STANDARD)
  string(REGEX MATCH [[..$]] VAL "${PYBIND11_CPP_STANDARD}")
  if(CMAKE_CXX_STANDARD)
    if(NOT CMAKE_CXX_STANDARD STREQUAL VAL)
      message(WARNING "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} does not match "
                      "PYBIND11_CPP_STANDARD=${PYBIND11_CPP_STANDARD}, "
                      "please remove PYBIND11_CPP_STANDARD from your cache")
    endif()
  else()
    set(supported_standards 11 14 17 20)
    if("${VAL}" IN_LIST supported_standards)
      message(WARNING "USE -DCMAKE_CXX_STANDARD=${VAL} instead of PYBIND11_CPP_STANDARD")
      set(CMAKE_CXX_STANDARD
          ${VAL}
          CACHE STRING "From PYBIND11_CPP_STANDARD")
    else()
      message(FATAL_ERROR "PYBIND11_CPP_STANDARD should be replaced with CMAKE_CXX_STANDARD "
                          "(last two chars: ${VAL} not understood as a valid CXX std)")
    endif()
  endif()
endif()

# --------------------- Python specifics -------------------------

# CMake 3.27 removes the classic FindPythonInterp if CMP0148 is NEW
if(CMAKE_VERSION VERSION_LESS "3.27")
  set(_pybind11_missing_old_python "OLD")
else()
  cmake_policy(GET CMP0148 _pybind11_missing_old_python)
endif()

# Check to see which Python mode we are in, new, old, or no python
if(PYBIND11_NOPYTHON)
  set(_pybind11_nopython ON)
  # We won't use new FindPython if PYBIND11_FINDPYTHON is defined and falselike
  # Otherwise, we use if FindPythonLibs is missing or if FindPython was already used
elseif(
  (NOT DEFINED PYBIND11_FINDPYTHON
   OR PYBIND11_FINDPYTHON STREQUAL "COMPAT"
   OR PYBIND11_FINDPYTHON)
  AND (_pybind11_missing_old_python STREQUAL "NEW"
       OR PYBIND11_FINDPYTHON STREQUAL "COMPAT"
       OR PYBIND11_FINDPYTHON
       OR Python_FOUND
       OR Python3_FOUND
      ))

  # New mode
  if(Python_FOUND OR Python3_FOUND)
    include("${CMAKE_CURRENT_LIST_DIR}/pybind11NewTools.cmake")
  else()
    include("${CMAKE_CURRENT_LIST_DIR}/pybind11NewTools.cmake")

    if(PYBIND11_FINDPYTHON STREQUAL "COMPAT")
      message(
        "Using compatibility mode for Python, set PYBIND11_FINDPYTHON to NEW/OLD to silence this message"
      )
      set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
      set(PYTHON_INCLUDE_DIR "${Python_INCLUDE_DIR}")
      set(Python_INCLUDE_DIRS "${Python_INCLUDE_DIRS}")
      set(PYTHON_LIBRARY "${Python_LIBRARY}")
      set(PYTHON_LIBRARIES "${Python_LIBRARIES}")
      set(PYTHON_VERSION "${Python_VERSION}")
      set(PYTHON_VERSION_STRING "${Python_VERSION_STRING}")
      set(PYTHON_VERSION_MAJOR "${Python_VERSION_MAJOR}")
      set(PYTHON_VERSION_MINOR "${Python_VERSION_MINOR}")
      set(PYTHON_VERSION_PATCH "${Python_VERSION_PATCH}")
    endif()
  endif()

else()

  # Classic mode
  include("${CMAKE_CURRENT_LIST_DIR}/pybind11Tools.cmake")

endif()

# --------------------- pybind11_find_import -------------------------------

if(NOT _pybind11_nopython AND NOT _PYBIND11_CROSSCOMPILING)
  # Check to see if modules are importable. Use REQUIRED to force an error if
  # one of the modules is not found. <package_name>_FOUND will be set if the
  # package was found (underscores replace dashes if present). QUIET will hide
  # the found message, and VERSION will require a minimum version. A successful
  # find will cache the result.
  function(pybind11_find_import PYPI_NAME)
    # CMake variables need underscores (PyPI doesn't care)
    string(REPLACE "-" "_" NORM_PYPI_NAME "${PYPI_NAME}")

    # Return if found previously
    if(${NORM_PYPI_NAME}_FOUND)
      return()
    endif()

    set(options "REQUIRED;QUIET")
    set(oneValueArgs "VERSION")
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "" ${ARGN})

    if(ARG_REQUIRED)
      set(status_level FATAL_ERROR)
    else()
      set(status_level WARNING)
    endif()

    execute_process(
      COMMAND
        ${${_Python}_EXECUTABLE} -c "
try:
    from importlib.metadata import version
except ImportError:
    from pkg_resources import get_distribution
    def version(s):
        return get_distribution(s).version
print(version('${PYPI_NAME}'))
        "
      RESULT_VARIABLE RESULT_PRESENT
      OUTPUT_VARIABLE PKG_VERSION
      ERROR_QUIET)

    string(STRIP "${PKG_VERSION}" PKG_VERSION)

    # If a result is present, this failed
    if(RESULT_PRESENT)
      set(${NORM_PYPI_NAME}_FOUND
          ${NORM_PYPI_NAME}-NOTFOUND
          CACHE INTERNAL "")
      # Always warn or error
      message(
        ${status_level}
        "Missing: ${PYPI_NAME} ${ARG_VERSION}\nTry: ${${_Python}_EXECUTABLE} -m pip install ${PYPI_NAME}"
      )
    else()
      if(ARG_VERSION AND PKG_VERSION VERSION_LESS ARG_VERSION)
        message(
          ${status_level}
          "Version incorrect: ${PYPI_NAME} ${PKG_VERSION} found, ${ARG_VERSION} required - try upgrading"
        )
      else()
        set(${NORM_PYPI_NAME}_FOUND
            YES
            CACHE INTERNAL "")
        set(${NORM_PYPI_NAME}_VERSION
            ${PKG_VERSION}
            CACHE INTERNAL "")
      endif()
      if(NOT ARG_QUIET)
        message(STATUS "Found ${PYPI_NAME} ${PKG_VERSION}")
      endif()
    endif()
    if(NOT ARG_VERSION OR (NOT PKG_VERSION VERSION_LESS ARG_VERSION))
      # We have successfully found a good version, cache to avoid calling again.
    endif()
  endfunction()
endif()

# --------------------- LTO -------------------------------

include(CheckCXXCompilerFlag)

# Checks whether the given CXX/linker flags can compile and link a cxx file.
# cxxflags and linkerflags are lists of flags to use.  The result variable is a
# unique variable name for each set of flags: the compilation result will be
# cached base on the result variable.  If the flags work, sets them in
# cxxflags_out/linkerflags_out internal cache variables (in addition to
# ${result}).
function(_pybind11_return_if_cxx_and_linker_flags_work result cxxflags linkerflags cxxflags_out
         linkerflags_out)
  set(CMAKE_REQUIRED_LIBRARIES ${linkerflags})
  check_cxx_compiler_flag("${cxxflags}" ${result})
  if(${result})
    set(${cxxflags_out}
        "${cxxflags}"
        PARENT_SCOPE)
    set(${linkerflags_out}
        "${linkerflags}"
        PARENT_SCOPE)
  endif()
endfunction()

function(_pybind11_generate_lto target prefer_thin_lto)
  if(MINGW)
    message(STATUS "${target} disabled (problems with undefined symbols for MinGW for now)")
    return()
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(cxx_append "")
    set(linker_append "")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT APPLE)
      # Clang Gold plugin does not support -Os; append -O3 to MinSizeRel builds to override it
      set(linker_append ";$<$<CONFIG:MinSizeRel>:-O3>")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND NOT MINGW)
      set(cxx_append ";-fno-fat-lto-objects")
    endif()

    if(prefer_thin_lto)
      set(thin "=thin")
    else()
      set(thin "")
    endif()

    if(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le" OR CMAKE_SYSTEM_PROCESSOR MATCHES "mips64")
      # Do nothing
    elseif(CMAKE_SYSTEM_NAME MATCHES Emscripten)
      # This compile is very costly when cross-compiling, so set this without checking
      set(PYBIND11_LTO_CXX_FLAGS "-flto${thin}${cxx_append}")
      set(PYBIND11_LTO_LINKER_FLAGS "-flto${thin}${linker_append}")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      _pybind11_return_if_cxx_and_linker_flags_work(
        HAS_FLTO_THIN "-flto${thin}${cxx_append}" "-flto${thin}${linker_append}"
        PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    endif()
    if(NOT HAS_FLTO_THIN)
      _pybind11_return_if_cxx_and_linker_flags_work(
        HAS_FLTO_AUTO "-flto=auto${cxx_append}" "-flto=auto${linker_append}"
        PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    endif()
    if(NOT HAS_FLTO_AUTO)
      _pybind11_return_if_cxx_and_linker_flags_work(
        HAS_FLTO "-flto${cxx_append}" "-flto${linker_append}" PYBIND11_LTO_CXX_FLAGS
        PYBIND11_LTO_LINKER_FLAGS)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    # IntelLLVM equivalent to LTO is called IPO; also IntelLLVM is WIN32/UNIX
    # WARNING/HELP WANTED: This block of code is currently not covered by pybind11 GitHub Actions!
    if(WIN32)
      _pybind11_return_if_cxx_and_linker_flags_work(
        HAS_INTEL_IPO "-Qipo" "-Qipo" PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    else()
      _pybind11_return_if_cxx_and_linker_flags_work(
        HAS_INTEL_IPO "-ipo" "-ipo" PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # Intel equivalent to LTO is called IPO
    _pybind11_return_if_cxx_and_linker_flags_work(HAS_INTEL_IPO "-ipo" "-ipo"
                                                  PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
  elseif(MSVC)
    # cmake only interprets libraries as linker flags when they start with a - (otherwise it
    # converts /LTCG to \LTCG as if it was a Windows path).  Luckily MSVC supports passing flags
    # with - instead of /, even if it is a bit non-standard:
    _pybind11_return_if_cxx_and_linker_flags_work(HAS_MSVC_GL_LTCG "/GL" "-LTCG"
                                                  PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
  endif()

  # Enable LTO flags if found, except for Debug builds
  if(PYBIND11_LTO_CXX_FLAGS)
    # CONFIG takes multiple values in CMake 3.19+, until then we have to use OR
    set(is_debug "$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>")
    set(not_debug "$<NOT:${is_debug}>")
    set(cxx_lang "$<COMPILE_LANGUAGE:CXX>")
    set(genex "$<AND:${not_debug},${cxx_lang}>")
    set_property(
      TARGET ${target}
      APPEND
      PROPERTY INTERFACE_COMPILE_OPTIONS "$<${genex}:${PYBIND11_LTO_CXX_FLAGS}>")
    if(CMAKE_PROJECT_NAME STREQUAL "pybind11")
      message(STATUS "${target} enabled")
    endif()
  else()
    if(CMAKE_PROJECT_NAME STREQUAL "pybind11")
      message(STATUS "${target} disabled (not supported by the compiler and/or linker)")
    endif()
  endif()

  if(PYBIND11_LTO_LINKER_FLAGS)
    set_property(
      TARGET ${target}
      APPEND
      PROPERTY INTERFACE_LINK_OPTIONS "$<${not_debug}:${PYBIND11_LTO_LINKER_FLAGS}>")
  endif()
endfunction()

if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  add_library(pybind11::lto IMPORTED INTERFACE ${optional_global})
  _pybind11_generate_lto(pybind11::lto FALSE)

  add_library(pybind11::thin_lto IMPORTED INTERFACE ${optional_global})
  _pybind11_generate_lto(pybind11::thin_lto TRUE)
endif()

# ---------------------- pybind11_strip -----------------------------

function(pybind11_strip target_name)
  # Strip unnecessary sections of the binary on Linux/macOS
  if(CMAKE_STRIP)
    if(APPLE)
      set(x_opt -x)
    endif()

    add_custom_command(
      TARGET ${target_name}
      POST_BUILD
      COMMAND ${CMAKE_STRIP} ${x_opt} $<TARGET_FILE:${target_name}>)
  endif()
endfunction()
