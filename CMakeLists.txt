cmake_minimum_required(VERSION 2.6)

project(LinearAlgebra)

set(EXECUTABLE_OUTPUT_PATH bin)

file(

        GLOB_RECURSE

        source_files

        src/*

)

add_executable(

        test
	
        ${source_files}

)
