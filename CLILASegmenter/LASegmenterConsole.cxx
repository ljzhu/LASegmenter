#include <iostream>
#include <stdlib.h>
#include "fstream"

//VTK
#include "vtkPolyDataReader.h"

// LASeg
//#include "./SegSrc/utilities.hxx"
#include "utilities.hxx"
#include "./SegSrc/LASeg.h"


int main(int argc, char** argv) {

    if(argc < 4) {
        std::cout << "missing arguments\n";
        return 0;
    }
    std::string sourceImageName(argv[1]);
    std::string seedImageName(argv[2]);
    std::string laImageName(argv[3]);

    LASeg::RSSWithShapeConstraint(sourceImageName.c_str(), seedImageName.c_str(), laImageName.c_str(), 0.1, 0.4);

    return EXIT_SUCCESS;

}
