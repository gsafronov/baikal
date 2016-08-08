set(PACKAGE breconstruct)

set(SRCS
  BReconstruct.cc
#  BReconstructLaser.cc
  BReconstructMuon.cc
  BRecParameters.cc
  BRecQualify.cc
  BRecQuality.cc
  )

include_directories(${ROOT_INCLUDE_DIR})
link_directories( ${ROOT_LIBRARY_DIR} )

String(REPLACE ".cc" ".h" HDRS "${SRCS}")
ROOT_GENERATE_DICTIONARY( G__${PACKAGE} ${HDRS} BRecIncl.h LINKDEF BRecLinkDef.h)
list(APPEND SRCS G__${PACKAGE}.cxx)

add_library (${PACKAGE} SHARED ${SRCS} ${HDRS} BRecIncl.h)
target_link_libraries (${PACKAGE} ${ROOT_LIBRARIES})

add_par(${PACKAGE})
