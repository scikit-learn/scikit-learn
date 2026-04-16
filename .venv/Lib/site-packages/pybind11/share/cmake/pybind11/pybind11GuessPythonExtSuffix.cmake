cmake_minimum_required(VERSION 3.15...4.2)

function(pybind11_guess_python_module_extension python)

  # The SETUPTOOLS_EXT_SUFFIX environment variable takes precedence:
  if(NOT DEFINED PYTHON_MODULE_EXT_SUFFIX AND DEFINED ENV{SETUPTOOLS_EXT_SUFFIX})
    message(
      STATUS
        "Getting Python extension suffix from ENV{SETUPTOOLS_EXT_SUFFIX}: $ENV{SETUPTOOLS_EXT_SUFFIX}"
    )
    set(PYTHON_MODULE_EXT_SUFFIX
        "$ENV{SETUPTOOLS_EXT_SUFFIX}"
        CACHE
          STRING
          "Extension suffix for Python extension modules (Initialized from SETUPTOOLS_EXT_SUFFIX)")
  endif()

  # The final extension depends on the system
  set(_PY_BUILD_EXTENSION "${CMAKE_SHARED_MODULE_SUFFIX}")
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(_PY_BUILD_EXTENSION ".pyd")
  endif()

  # If running under scikit-build-core, use the SKBUILD_SOABI variable:
  if(NOT DEFINED PYTHON_MODULE_EXT_SUFFIX AND DEFINED SKBUILD_SOABI)
    message(STATUS "Determining Python extension suffix based on SKBUILD_SOABI: ${SKBUILD_SOABI}")
    set(PYTHON_MODULE_EXT_SUFFIX ".${SKBUILD_SOABI}${_PY_BUILD_EXTENSION}")
  endif()

  # If that didn't work, use the Python_SOABI variable:
  if(NOT DEFINED PYTHON_MODULE_EXT_SUFFIX AND DEFINED ${python}_SOABI)
    message(
      STATUS "Determining Python extension suffix based on ${python}_SOABI: ${${python}_SOABI}")
    # If the SOABI already has an extension, use it as the full suffix
    # (used for debug versions of Python on Windows)
    if(${python}_SOABI MATCHES "\\.")
      set(PYTHON_MODULE_EXT_SUFFIX "${${python}_SOABI}")
      # If the SOABI is empty, this is usually a bug, but we generate a
      # correct extension anyway, which is the best we can do
    elseif("${${python}_SOABI}" STREQUAL "")
      message(
        WARNING
          "${python}_SOABI is defined but empty. You may want to set PYTHON_MODULE_EXT_SUFFIX explicitly."
      )
      set(PYTHON_MODULE_EXT_SUFFIX "${_PY_BUILD_EXTENSION}")
      # Otherwise, add the system-dependent extension to it
    else()
      set(PYTHON_MODULE_EXT_SUFFIX ".${${python}_SOABI}${_PY_BUILD_EXTENSION}")
    endif()
  endif()

  # If we could not deduce the extension suffix, unset the results:
  if(NOT DEFINED PYTHON_MODULE_EXT_SUFFIX)
    unset(PYTHON_MODULE_DEBUG_POSTFIX CACHE)
    unset(PYTHON_MODULE_EXTENSION CACHE)
    unset(PYTHON_IS_DEBUG CACHE)
    return()
  endif()

  # Sanity checks:
  if(${python}_SOABI AND NOT (PYTHON_MODULE_EXT_SUFFIX STREQUAL ${python}_SOABI
                              OR PYTHON_MODULE_EXT_SUFFIX MATCHES "\\.${${python}_SOABI}\\."))
    message(
      WARNING
        "Python extension suffix (${PYTHON_MODULE_EXT_SUFFIX}) does not match ${python}_SOABI (${${python}_SOABI})."
    )
  endif()

  # Separate file name postfix from extension: (https://github.com/pybind/pybind11/issues/4699)
  get_filename_component(_PYTHON_MODULE_DEBUG_POSTFIX "${PYTHON_MODULE_EXT_SUFFIX}" NAME_WE)
  get_filename_component(_PYTHON_MODULE_EXTENSION "${PYTHON_MODULE_EXT_SUFFIX}" EXT)

  # Try to deduce the debug ABI from the extension suffix:
  if(NOT DEFINED _PYTHON_IS_DEBUG)
    if(_PYTHON_MODULE_EXTENSION MATCHES "^\\.(cpython-|cp|pypy)[0-9]+dm?-"
       OR _PYTHON_MODULE_DEBUG_POSTFIX MATCHES "^_d")
      set(_PYTHON_IS_DEBUG On)
    else()
      set(_PYTHON_IS_DEBUG Off)
    endif()
  endif()

  # Return results
  set(PYTHON_MODULE_DEBUG_POSTFIX
      "${_PYTHON_MODULE_DEBUG_POSTFIX}"
      CACHE INTERNAL "")
  set(PYTHON_MODULE_EXTENSION
      "${_PYTHON_MODULE_EXTENSION}"
      CACHE INTERNAL "")
  set(PYTHON_IS_DEBUG
      "${_PYTHON_IS_DEBUG}"
      CACHE INTERNAL "")

endfunction()
