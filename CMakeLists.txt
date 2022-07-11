cmake_minimum_required( VERSION 3.5 FATAL_ERROR )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3" )

message( "Compiler: ${CMAKE_CXX_COMPILER}" )
message( "Option: ${CMAKE_CXX_FLAGS}" )

# preprocessor macros
add_compile_definitions( NDEBUG )

include_directories( "." )

# phat
include_directories( "../phat-master/include/" )

# Boost
find_package( Boost 1.5 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} )

add_executable( dpc_vine
    dpc_main.cpp
    total_complex.cpp
    utils.cpp
    common.cpp
    dpc.cpp
    chain.cpp
    ud.cpp
    dynamic_zigzag.cpp
    fzz.cpp
    ) 

target_link_libraries( dpc_vine 
    ${Boost_PROGRAM_OPTIONS_LIBRARY} 
    )
