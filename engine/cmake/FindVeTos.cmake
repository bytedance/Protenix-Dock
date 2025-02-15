# Cache Variables: (probably not for direct use in your scripts)
#  VeTos_INCLUDE_DIR
#  VeTos_LIBRARY
#
# Non-cache variables you might use in your CMakeLists.txt:
#  VeTos_FOUND
#  VeTos_INCLUDE_DIRS
#  VeTos_LIBRARIES
#
# Requires these CMake modules:
#  FindPackageHandleStandardArgs (known included with CMake >=2.6.2)

set(VeTos_ROOT_DIR "/opt/tiger/ve-tos-cpp-sdk" CACHE PATH "Root directory of VeTosCppSdk files.")

find_library(VeTos_LIBRARY NAMES ve-tos-cpp-sdk-lib
                           PATHS "${VeTos_ROOT_DIR}"
                           PATH_SUFFIXES "lib/")
get_filename_component(_libdir "${VeTos_LIBRARY}" DIRECTORY)
find_path(VeTos_INCLUDE_DIR NAMES TosClientV2.h
                            HINTS "${_libdir}/.."
                            PATHS "${VeTos_ROOT_DIR}"
                            PATH_SUFFIXES "include/")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VeTos DEFAULT_MSG VeTos_LIBRARY VeTos_INCLUDE_DIR)

if(VeTos_FOUND)
    set(VeTos_LIBRARIES "${VeTos_LIBRARY}")
    set(VeTos_INCLUDE_DIRS "${VeTos_INCLUDE_DIR}")
    mark_as_advanced(VeTos_ROOT_DIR)
endif()
mark_as_advanced(VeTos_LIBRARY VeTos_INCLUDE_DIR)
