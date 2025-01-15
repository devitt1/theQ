cmake_minimum_required(VERSION 3.10)
project(StimTest)

set(STIM_DIR "./Stim")


include_directories(${STIM_DIR})


set(STIM_LIB_PATH "${STIM_DIR}/build/out")


add_executable(stim_test stim_test.cpp)


target_link_libraries(stim_test ${STIM_LIB_PATH}/libstim.a)
