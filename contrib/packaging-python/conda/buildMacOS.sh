
mkdir -p ./build
cd ./build
if [ "$PY3K" == "1" ]; then
    MY_PY_VER="${PY_VER}m"
else
    MY_PY_VER="${PY_VER}"
fi

if [ `uname` == Darwin ]; then
    PY_LIB="libpython${MY_PY_VER}.dylib"
else
    PY_LIB="libpython${MY_PY_VER}.so"
fi

# set MKL vars
export MKL_INTERFACE_LAYER=LP64
export MKL_THREADING_LAYER=INTEL

if [ `uname` == Darwin ]; then
    sed -i '' 's/${PYTHON_LIBRARY}//g' $SRC_DIR/src/chrono_python/CMakeLists.txt
    sed -i '' 's/find_package(AVX)//g' $SRC_DIR/src/CMakeLists.txt
    sed -i '' 's/find_package(SSE)//g' $SRC_DIR/src/CMakeLists.txt
    sed -i '' 's/find_package(NEON)//g' $SRC_DIR/src/CMakeLists.txt
    sed -i '' 's/find_package(FMA)//g' $SRC_DIR/src/CMakeLists.txt
    sed -i '' 's/find_package(OpenMP)//g' $SRC_DIR/src/CMakeLists.txt
fi
export LDFLAGS="-Wl,-undefined,dynamic_lookup $LDFLAGS"

CONFIGURATION=Release
# Configure step
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
 -DCMAKE_C_COMPILER=$(which clang) \
 -DCMAKE_CXX_COMPILER=$(which clang++) \
 -DCMAKE_PREFIX_PATH=$PREFIX \
 -DCMAKE_SYSTEM_PREFIX_PATH=$PREFIX \
 -DCH_INSTALL_PYTHON_PACKAGE=$SP_DIR \
 -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON \
 -DPYTHON_INCLUDE_DIR:PATH=$PREFIX/include/python$MY_PY_VER \
 -DPYTHON_LIBRARY:FILEPATH=$PREFIX/lib/${PY_LIB} \
 -DCMAKE_BUILD_TYPE=$CONFIGURATION \
 -DENABLE_MODULE_IRRLICHT=ON \
 -DENABLE_MODULE_POSTPROCESS=ON \
 -DENABLE_MODULE_VEHICLE=ON \
 -DENABLE_MODULE_PYTHON=ON \
 -DBUILD_DEMOS=OFF \
 -DBUILD_TESTING=OFF \
 -DBUILD_BENCHMARKING=OFF \
 -DBUILD_GMOCK=OFF \
 -DENABLE_MODULE_CASCADE=OFF \
 -DCASCADE_INCLUDE_DIR=$HOME/miniconda/include/opencascade \
 -DCASCADE_LIBDIR=$HOME/miniconda/lib \
 -DENABLE_MODULE_MKL=OFF \
 -DMKL_INCLUDE_DIR=$HOME/miniconda/include \
 -DMKL_RT_LIBRARY=$HOME/miniconda/lib/libmkl_rt.a \
 -DEIGEN3_INCLUDE_DIR=/usr/local/include/eigen3 \
 -DPYCHRONO_DATA_PATH=../../../../../../share/chrono/data \
 ./..
# Build step
# on linux travis, limit the number of concurrent jobs otherwise
# gcc gets out of memory
cmake --build . --config "$CONFIGURATION"

cmake --build . --config "$CONFIGURATION" --target install



