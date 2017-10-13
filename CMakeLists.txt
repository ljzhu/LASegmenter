cmake_minimum_required(VERSION 3.5)

project(LASegmenter)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php?title=Documentation/Nightly/Extensions/LASegmenter")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Liangjia Zhu (SBU), Yi Gao (UAB), Josh Cates (Utah), Alan Morris (Utah), Danny Perry (Utah), Greg Gardner (Utah), Rob MacLeod (Utah), Allen Tannenbaum (SBU) ")
set(EXTENSION_DESCRIPTION "This module performs segmentation of the left atrium from MR images")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/0/0b/LASegmenter.png")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/a/a7/LeftAtriumSegmenterScreenShot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CLILASegmenter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
