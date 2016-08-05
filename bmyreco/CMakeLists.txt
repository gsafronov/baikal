set(PACKAGE bmyreco)

set(SRCS
BMyReco.cc
BMyRecoMC.cc
BMyRecoBEvt.cc
BMyRecoReco.cc

include_directories(${ROOT_INCLUDE_DIR})
link_directories( ${ROOT_LIBRARY_DIR} )

String(REPLACE ".cc" ".h" HDRS "${SRCS}")
ROOT_GENERATE_DICTIONARY( G__${PACKAGE} ${HDRS} LINKDEF BMyRecLinkDef.h)
list(APPEND SRCS G__${PACKAGE}.cxx)

add_library (${PACKAGE} SHARED ${SRCS} ${HDRS})
target_link_libraries (${PACKAGE} ${ROOT_LIBRARIES})

add_par(${PACKAGE})
