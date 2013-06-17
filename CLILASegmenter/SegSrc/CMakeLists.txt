set(LASegSrc_Files LASeg.txx)
add_library(SegSrc ${LASegSrc_Files})
set_target_properties(SegSrc PROPERTIES LINKER_LANGUAGE CXX)
