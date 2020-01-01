set( analyze_SOURCE
  anal_qn.cc
  anal_PHENIXspectra.cc
  anal_STARspectra.cc
  anal_balance.cc
  anal_gamma.cc
  anal_hbt.cc
  anal_v2.cc
  analyze.cc
  consolidate.cc
  distill.cc
)

include_directories(
  ${coral_INCLUDE_DIRS}

  ../../run/
)

add_library( analyze ${analyze_SOURCE} )

install( TARGETS analyze
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install( FILES
    analyze.h
  DESTINATION
    include
)
