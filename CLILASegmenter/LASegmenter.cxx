#include "LASegmenterCLP.h"
#include <iostream>
#include "fstream"

//VTK
#include "vtkPolyDataReader.h"

// LASeg
#include "utilities.hxx"
#include "./SegSrc/LASeg.h"

int main(int argc, char** argv) {

     PARSE_ARGS;

     LASeg::RSSWithShapeConstraint(sourceImageName.c_str(), seedImageName.c_str(), laImageName.c_str(), intensityHomogeneity, curvatureWeight);

    return EXIT_SUCCESS;
}
