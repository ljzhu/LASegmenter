set(ShapeMoments_Src Binomial.h Binomial.cpp Factorial.h ScaledGeometricMoments.cpp ZernikeDescriptor.cpp ZernikeMoments.cpp)
add_library(ShapeMoments ${ShapeMoments_Src})
set_target_properties(ShapeMoments PROPERTIES LINKER_LANGUAGE CXX)
