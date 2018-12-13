IF(COPY_PARAMFILES)
  SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/ConvectionDiffusion2D)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.gup ${_SRC}/*.bup .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE("ConvectionDiffusion2D${BUILD_ENDING}" main.cc local.cc localequation.cc)

