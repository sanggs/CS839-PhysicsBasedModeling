add_executable (assignment2
        main.cpp
        FiniteElementMesh.h
        ConjugateGradient.h
        CGSystemWrapper.h
        CGVectorWrapper.h
        AnimatedTetrahedronMesh.h
        )

find_package(PythonLibs)

set(EIGEN3_INC_DIR /usr/local/include/eigen3) # Eigen header search path in mac (default). Change this if different in your system.

target_include_directories(assignment2
        PUBLIC
            ${USD_INC_DIR}
            ${EIGEN3_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
        )

target_link_libraries(
        assignment2
                ${USD_LIB_DIR}/libgf.dylib
                ${USD_LIB_DIR}/libsdf.dylib
                ${USD_LIB_DIR}/libtf.dylib
                ${USD_LIB_DIR}/libusd.dylib
                ${USD_LIB_DIR}/libusdGeom.dylib
                ${USD_LIB_DIR}/libvt.dylib
                ${USD_LIB_DIR}/libboost_python.dylib # todo find library
                ${PYTHON_LIBRARY}
        )
