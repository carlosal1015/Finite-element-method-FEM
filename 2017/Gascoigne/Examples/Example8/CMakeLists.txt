IF(COPY_PARAMFILES)
  SET(_SRC ${GASCOIGNESTD_SOURCE_DIR}/Examples/Example8)
  ADD_CUSTOM_TARGET(COPY ALL cp ${_SRC}/*.param ${_SRC}/*.inp .)
ENDIF(COPY_PARAMFILES)

ADD_EXECUTABLE("Gascoigne8${BUILD_ENDING}" navierstokessplit.cc chorinalgorithm.cc main.cc)



