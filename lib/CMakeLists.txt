add_library(libs STATIC initialize.cpp hamiltonian.cpp fermi_surface.cpp band.cpp green_function.cpp response.cpp)
target_include_directories(libs PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(libs PRIVATE CGAL::CGAL)
