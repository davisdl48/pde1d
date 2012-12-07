project(pde1d)
cmake_minimum_required(VERSION 2.6)
find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR})

set(pde1d_SRCS pde1d.cpp main.cpp)
qt4_automoc(${pde1d_SRCS})
add_executable(pde1d ${pde1d_SRCS})
target_link_libraries(pde1d ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})
