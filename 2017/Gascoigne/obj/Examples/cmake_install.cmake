# Install script for directory: /home/carlosal1015/Gascoigne/Examples

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
    set(CMAKE_INSTALL_CONFIG_NAME "")
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
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example1/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example2/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example3/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example4/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example5/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example6/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example7/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example8/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Example9/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/MeshTest/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Laplace2D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Laplace3D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/NavierStokes2D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/NavierStokes3D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/NavierStokesFace/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Heat2D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/ConvectionDiffusion2D/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/Functionals/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/MeshInterpolator/cmake_install.cmake")
  include("/home/carlosal1015/Gascoigne/obj/Examples/PeriodicStokes2D/cmake_install.cmake")

endif()

