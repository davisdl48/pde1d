CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#FIND_PACKAGE(Qt4 REQUIRED)

INCLUDE_DIRECTORIES(/usr/include/qwt6 ${CMAKE_BINARY_DIR}/src ${CMAKE_BINARY_DIR}/superlu)

INCLUDE (  ${CMAKE_ROOT}/Modules/FindQt4.cmake  )
INCLUDE(${QT_USE_FILE})

FILE( GLOB pde1d_RESOURCES *.qrc)
QT4_ADD_RESOURCES(pde1d_RESOURCES_RCC ${pde1d_RESOURCES})

FILE( GLOB pde1d_FORMS *.ui)
QT4_WRAP_UI(pde1d_FORMS_HEADERS ${pde1d_FORMS})

set( SOLVS pswidget.h rkwidget.h envwidget.h specwidget.h impwidget.h femwidget.h  leastsqrwidget.h simpimpwidget.h idealwidget.h)
set( pde1d_MOC_HEADERS pde1d.h solvwidget.h controls.h errtabdock.h myinputs.h curvetabdock.h itemdelegate.h curvesmodel.h colordelegate.h ${SOLVS} )
# 
QT4_WRAP_CPP(pde1d_HEADERS_MOC ${pde1d_MOC_HEADERS})

FILE( GLOB pde1d_SRCS *.cpp)
FILE( GLOB pde1d_HEADERS *.h)
#QT4_AUTOMOC(${pde1d_SRCS} ${pde1d_MOC_HEADERS})

ADD_EXECUTABLE(pde1d ${pde1d_SRCS} ${pde1d_HEADERS_MOC} ${pde1d_FORMS_HEADERS} ${pde1d_RESOURCES_RCC}  )
# ../superlu/lib/libsuperlu_4.3.a)
TARGET_LINK_LIBRARIES(pde1d ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${CMAKE_BINARY_DIR}/superlu/libsuperlu.a gsl gslcblas m qwt blas)

option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
      "Doxygen is needed to build the documentation. Please install it correctly")
  endif()
  #-- Configure the Template Doxyfile for our specific project
  configure_file(Doxyfile.in 
                 ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
  #-- Add a custom target to run Doxygen when ever the project is built
  add_custom_target (Docs ALL 
  COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
  SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
  # IF you do NOT want the documentation to be generated EVERY time you build the project
  # then leave out the 'ALL' keyword from the above command.
endif()
