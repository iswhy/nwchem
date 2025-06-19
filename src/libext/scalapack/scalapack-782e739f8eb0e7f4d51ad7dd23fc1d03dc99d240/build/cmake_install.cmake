# Install script for directory: /home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "MinSizeRel")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/lib/libscalapack.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/pkgconfig/scalapack.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/lib/pkgconfig" TYPE FILE FILES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/scalapack.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1" TYPE FILE FILES
    "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CMakeFiles/scalapack-config.cmake"
    "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/scalapack-config-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1/scalapack-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1/scalapack-targets.cmake"
         "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CMakeFiles/Export/51f8d17d98dae6d66dec5275f1680311/scalapack-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1/scalapack-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1/scalapack-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1" TYPE FILE FILES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CMakeFiles/Export/51f8d17d98dae6d66dec5275f1680311/scalapack-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/scalapack-2.2.1" TYPE FILE FILES "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/CMakeFiles/Export/51f8d17d98dae6d66dec5275f1680311/scalapack-targets-minsizerel.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/BLACS/cmake_install.cmake")
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/TOOLS/cmake_install.cmake")
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/PBLAS/cmake_install.cmake")
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/REDIST/cmake_install.cmake")
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/SRC/cmake_install.cmake")
  include("/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/TESTING/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/imam77/src/nwchem/NWChem-WFAT/src/libext/scalapack/scalapack/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
