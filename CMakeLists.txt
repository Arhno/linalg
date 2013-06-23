cmake_minimum_required(VERSION 2.6)

project(LinearAlgebra)

file(

        GLOB_RECURSE

        header_files

        include/*

)

file(

        GLOB_RECURSE

        source_files

        src/*

)

add_executable(

        test
	
	${header_files}
        ${source_files}

)
