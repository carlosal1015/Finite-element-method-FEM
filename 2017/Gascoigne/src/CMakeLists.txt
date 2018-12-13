AUX_SOURCE_DIRECTORY(Common Common_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Q1 DiscretizationQ1_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Q2 DiscretizationQ2_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Q4 DiscretizationQ4_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Gls DiscretizationGls_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Lps DiscretizationLps_SOURCES)
AUX_SOURCE_DIRECTORY(Discretization/Face DiscretizationFace_SOURCES)
AUX_SOURCE_DIRECTORY(Solver Solver_SOURCES)
AUX_SOURCE_DIRECTORY(Problems Problems_SOURCES)
AUX_SOURCE_DIRECTORY(Mesh Mesh_SOURCES)
AUX_SOURCE_DIRECTORY(LinAlg LinAlg_SOURCES)

SET(GASCOIGNE_SOURCES 
${Common_SOURCES} 
${DiscretizationQ1_SOURCES} 
${DiscretizationQ2_SOURCES} 
${DiscretizationQ4_SOURCES} 
${DiscretizationGls_SOURCES} 
${DiscretizationLps_SOURCES} 
${DiscretizationFace_SOURCES} 
${Solver_SOURCES} 
${Problems_SOURCES} 
${LinAlg_SOURCES} 
${Mesh_SOURCES}
)

ADD_LIBRARY("GascoigneStd${BUILD_ENDING}" ${GASCOIGNE_SOURCES})

IF(APPLE)
	TARGET_LINK_LIBRARIES("GascoigneStd${BUILD_ENDING}" "-framework vecLib")
ELSE(APPLE)
	SET(GASCOIGNE_LIBS "GascoigneStd${BUILD_ENDING}" metis rt)
ENDIF(APPLE)
