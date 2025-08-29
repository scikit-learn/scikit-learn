# fail noisily if attempt to use this file without setting:
# cmake_minimum_required(VERSION ${CMAKE_VERSION})
# project(... LANGUAGES ...)

cmake_policy(SET CMP0000 NEW)

set(PACKAGE_FOUND FALSE)
set(_packageName "${NAME}")
string(TOUPPER "${_packageName}" PACKAGE_NAME)

if("${STATIC}" STREQUAL "True")
  set("${NAME}_USE_STATIC_LIBS" "ON")
endif()

while(TRUE)
  if ("${VERSION}" STREQUAL "")
    find_package("${NAME}" QUIET COMPONENTS ${COMPS})
  else()
    find_package("${NAME}" "${VERSION}" QUIET COMPONENTS ${COMPS})
  endif()

  # ARCHS has to be set via the CMD interface
  if(${_packageName}_FOUND OR ${PACKAGE_NAME}_FOUND OR "${ARCHS}" STREQUAL "")
    break()
  endif()

  list(GET       ARCHS 0 CMAKE_LIBRARY_ARCHITECTURE)
  list(REMOVE_AT ARCHS 0)
endwhile()

if(${_packageName}_FOUND  OR  ${PACKAGE_NAME}_FOUND)
  set(PACKAGE_FOUND TRUE)

  # Check the following variables:
  # FOO_VERSION
  # Foo_VERSION
  # FOO_VERSION_STRING
  # Foo_VERSION_STRING
  if(NOT DEFINED PACKAGE_VERSION)
    if(DEFINED ${_packageName}_VERSION)
      set(PACKAGE_VERSION "${${_packageName}_VERSION}")
    elseif(DEFINED ${PACKAGE_NAME}_VERSION)
      set(PACKAGE_VERSION "${${PACKAGE_NAME}_VERSION}")
    elseif(DEFINED ${_packageName}_VERSION_STRING)
      set(PACKAGE_VERSION "${${_packageName}_VERSION_STRING}")
    elseif(DEFINED ${PACKAGE_NAME}_VERSION_STRING)
      set(PACKAGE_VERSION "${${PACKAGE_NAME}_VERSION_STRING}")
    endif()
  endif()

  # Check the following variables:
  # FOO_LIBRARIES
  # Foo_LIBRARIES
  # FOO_LIBS
  # Foo_LIBS
  set(libs)
  if(DEFINED ${_packageName}_LIBRARIES)
    set(libs ${_packageName}_LIBRARIES)
  elseif(DEFINED ${PACKAGE_NAME}_LIBRARIES)
    set(libs ${PACKAGE_NAME}_LIBRARIES)
  elseif(DEFINED ${_packageName}_LIBS)
    set(libs ${_packageName}_LIBS)
  elseif(DEFINED ${PACKAGE_NAME}_LIBS)
    set(libs ${PACKAGE_NAME}_LIBS)
  endif()

  # Check the following variables:
  # FOO_INCLUDE_DIRS
  # Foo_INCLUDE_DIRS
  # FOO_INCLUDES
  # Foo_INCLUDES
  # FOO_INCLUDE_DIR
  # Foo_INCLUDE_DIR
  set(includes)
  if(DEFINED ${_packageName}_INCLUDE_DIRS)
    set(includes ${_packageName}_INCLUDE_DIRS)
  elseif(DEFINED ${PACKAGE_NAME}_INCLUDE_DIRS)
    set(includes ${PACKAGE_NAME}_INCLUDE_DIRS)
  elseif(DEFINED ${_packageName}_INCLUDES)
    set(includes ${_packageName}_INCLUDES)
  elseif(DEFINED ${PACKAGE_NAME}_INCLUDES)
    set(includes ${PACKAGE_NAME}_INCLUDES)
  elseif(DEFINED ${_packageName}_INCLUDE_DIR)
    set(includes ${_packageName}_INCLUDE_DIR)
  elseif(DEFINED ${PACKAGE_NAME}_INCLUDE_DIR)
    set(includes ${PACKAGE_NAME}_INCLUDE_DIR)
  endif()

  # Check the following variables:
  # FOO_DEFINITIONS
  # Foo_DEFINITIONS
  set(definitions)
  if(DEFINED ${_packageName}_DEFINITIONS)
    set(definitions ${_packageName}_DEFINITIONS)
  elseif(DEFINED ${PACKAGE_NAME}_DEFINITIONS)
    set(definitions ${PACKAGE_NAME}_DEFINITIONS)
  endif()

  set(PACKAGE_INCLUDE_DIRS "${${includes}}")
  set(PACKAGE_DEFINITIONS  "${${definitions}}")
  set(PACKAGE_LIBRARIES    "${${libs}}")
endif()

if(${_packageName}_NOT_FOUND_MESSAGE)
  set(PACKAGE_NOT_FOUND_MESSAGE "${${_packageName}_NOT_FOUND_MESSAGE}")
elseif(${PACKAGE_NAME}_NOT_FOUND_MESSAGE)
  set(PACKAGE_NOT_FOUND_MESSAGE "${${PACKAGE_NAME}_NOT_FOUND_MESSAGE}")
endif()
