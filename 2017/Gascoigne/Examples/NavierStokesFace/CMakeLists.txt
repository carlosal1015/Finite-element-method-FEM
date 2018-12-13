IF(COPY_PARAMFILES)
  SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/NavierStokesFace)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE("NavierStokesFace${BUILD_ENDING}" main.cc nsface.cc)


