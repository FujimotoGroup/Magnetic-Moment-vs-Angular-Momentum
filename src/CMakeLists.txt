add_executable(main main.cpp)
target_link_libraries(main PUBLIC libs)
target_link_libraries(main PUBLIC mkl_intel_lp64)
target_link_libraries(main PUBLIC mkl_sequential)
target_link_libraries(main PUBLIC mkl_core)
target_link_libraries(main PUBLIC pthread)

target_compile_options(main PUBLIC -O1 -frounding-math)


