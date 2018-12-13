FIND_PACKAGE(Doxygen)
if (DOXYGEN_FOUND)
  # Set up tag file generation and store location in source file for other projects (like RoDoBo)
  set (DOXYGEN_TAG_FILE_OUTPUT "${PROJECT_BINARY_DIR}/doc/${PROJECT_NAME}.tag")
  # Create a CMake file containing the location of the tag file 
  configure_file(taginformation.cmake.in ${PROJECT_SOURCE_DIR}/doc/taginformation.cmake @ONLY)

  # set up additional tag files (except for libstdc++) to include
  set (DOXYGEN_ADDITIONAL_TAGFILES "")
  # include ("${GASCOIGNE_STD}/doc/taginformation.cmake" OPTIONAL)

  # Configure the Doxyfile
  configure_file(Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/doc/Doxyfile  @ONLY)
  # Add a custom target doc (has to be called by hand)
  add_custom_target (doc 
                    COMMAND ${DOXYGEN_EXECUTABLE} 
                    ${PROJECT_BINARY_DIR}/doc/Doxyfile 
                    SOURCES ${PROJECT_BINARY_DIR}/doc/Doxyfile
                     )
else(DOXYGEN_FOUND)
   message(STATUS "Doxygen not found, make doc not available.")
endif()
