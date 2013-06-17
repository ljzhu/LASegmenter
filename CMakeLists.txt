cmake_minimum_required(VERSION 2.8.9)

project(LASegmenter)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/LASegmenter")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "LiangJia Zhu (UAB), Yi Gao (BWH), Josh Cates (Utah), Alan Morris (Utah), Danny Perry (Utah), Greg Gardner (Utah), Rob MacLeod (Utah), Allen Tannenbaum (UAB) ")
set(EXTENSION_DESCRIPTION "This module performs segmentation of the left atrium from MR images")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/LASegmenter/LASegmenter.png?revision=21745&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-LASegmenter-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-LASegmenter-screenshot-2.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CLILASegmenter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
