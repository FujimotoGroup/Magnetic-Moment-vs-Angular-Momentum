add_executable(main main.cpp)
target_link_libraries(main libs)
target_link_libraries(main mkl_intel_lp64)
target_link_libraries(main mkl_sequential)
target_link_libraries(main mkl_core)

target_compile_options(main PUBLIC -O2 -Wall)


