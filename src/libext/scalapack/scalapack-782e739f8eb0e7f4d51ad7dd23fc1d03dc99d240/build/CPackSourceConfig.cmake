# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BUILD_SOURCE_DIRS "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack;/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/libext_utils/cmake-3.24.0-linux-x86_64/share/cmake-3.24/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "SCALAPACK built using CMake")
set(CPACK_DMG_SLA_USE_RESOURCE_FILE_LICENSE "ON")
set(CPACK_GENERATOR "TGZ")
set(CPACK_IGNORE_FILES "~$;.svn")
set(CPACK_INSTALLED_DIRECTORIES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack;/")
set(CPACK_INSTALL_CMAKE_PROJECTS "")
set(CPACK_INSTALL_PREFIX "/usr/local")
set(CPACK_MODULE_PATH "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/CMAKE")
set(CPACK_NSIS_DISPLAY_NAME "SCALAPACK")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_PACKAGE_NAME "SCALAPACK")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_OUTPUT_CONFIG_FILE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/libext_utils/cmake-3.24.0-linux-x86_64/share/cmake-3.24/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ScaLAPACK- Linear Algebra Package")
set(CPACK_PACKAGE_FILE_NAME "scalapack-2.2.1")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "SCALAPACK")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "SCALAPACK")
set(CPACK_PACKAGE_NAME "ScaLAPACK")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "University of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd")
set(CPACK_PACKAGE_VERSION "2.2.1")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/LICENSE")
set(CPACK_RESOURCE_FILE_README "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/libext_utils/cmake-3.24.0-linux-x86_64/share/cmake-3.24/Templates/CPack.GenericDescription.txt")
set(CPACK_RESOURCE_FILE_WELCOME "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/libext_utils/cmake-3.24.0-linux-x86_64/share/cmake-3.24/Templates/CPack.GenericWelcome.txt")
set(CPACK_RPM_PACKAGE_SOURCES "ON")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "~$;.svn")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack;/")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "scalapack-2.2.1")
set(CPACK_SOURCE_TOPLEVEL_TAG "Linux-Source")
set(CPACK_STRIP_FILES "")
set(CPACK_SYSTEM_NAME "Linux")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "Linux-Source")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
