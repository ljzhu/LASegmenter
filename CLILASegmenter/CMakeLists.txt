
cmake_minimum_required(VERSION 3.5)

SEMMacroBuildCLI(
  NAME LASegmenter
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${PROJECT_SOURCE_DIR}/CLILASegmenter/
    ${PROJECT_SOURCE_DIR}/CLILASegmenter/ShapeMoments
    ${PROJECT_SOURCE_DIR}/CLILASegmenter/lib
    ${PROJECT_SOURCE_DIR}/CLILASegmenter/RSSZM
    ${PROJECT_SOURCE_DIR}/CLILASegmenter/SegSrc
    ${PROJECT_BINARY_DIR} # add the binary tree to the search path for include files
                          # so that we will find TutorialConfig.h
  #EXECUTABLE_ONLY
  )

