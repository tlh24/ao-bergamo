cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(torch_dmctrl)

find_package(Torch REQUIRED)

include_directories("/usr/local/include" "../../include/")# for cnpy
link_directories("/usr/local/lib" "../../AlpaoLib/Linux/x64/")

add_executable(torch_dmctrl main.cpp gettime.cpp sock.cpp)
target_link_libraries(torch_dmctrl cnpy z asdk gsl "${TORCH_LIBRARIES}" )
set_property(TARGET torch_dmctrl PROPERTY CXX_STANDARD 14)

# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_VERBOSE_MAKEFILE ON)

# after manually installing CUDNN, 
# run this with : 
# cmake src/ -DCMAKE_PREFIX_PATH="/home/tlh24/Dropbox/var/libtorch/"

# to clean: 
# rm -rf CMake* cmake_* detect_* Makefile

# prerequisites: 
# sudo apt-get install cmake
# download cudnn dev from nvidia (deb package works fine)
# install cnpy from git

