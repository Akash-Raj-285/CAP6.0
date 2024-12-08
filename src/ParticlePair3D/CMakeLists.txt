################################################################################################
# Project CAP/ParticlePair3D
################################################################################################
#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

ROOT_GENERATE_DICTIONARY(G__ParticlePair3D ParticlePair3DHistos.hpp ParticlePair3DDerivedHistos.hpp ParticlePair3DAnalyzer.hpp ParticlePair3DCalculator.hpp ParticlePair3DBfCalculator.hpp ParticlePair3DBfHistos.hpp LINKDEF ParticlePair3DLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -pedantic)
add_library(ParticlePair3D SHARED ParticlePair3DHistos.cpp ParticlePair3DDerivedHistos.cpp  ParticlePair3DAnalyzer.cpp  ParticlePair3DCalculator.cpp ParticlePair3DBfCalculator.cpp ParticlePair3DBfHistos.cpp G__ParticlePair3D.cxx)

target_link_libraries(ParticlePair3D Base ParticleDb  Particles ParticleSingle ParticlePair ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(ParticlePair3D  PUBLIC Base Particles ParticleSingle ParticlePair ParticlePair3D ${EXTRA_INCLUDES} )

# optimization for big histogram access within the pair inner loop
# WARNING: the sumw2 structure is not used -> estimate uncertainties using sub-samples
add_definitions(-DOPTIMIZEADDBINCONTENT)


install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libParticlePair3D.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libParticlePair3D_rdict.pcm" DESTINATION "$ENV{CAP_LIB_PATH}")
install(TARGETS ParticlePair3D  LIBRARY DESTINATION "$ENV{CAP_LIB_PATH}")

