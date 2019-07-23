file(GLOB NN_INCLUDE_NETWORKS "include/Networks/*.hpp")
file(GLOB NN_INCLUDE_OPTIMIZERS "include/Optimizers/*.hpp")
file(GLOB NN_SRC "src/*.cpp")

source_group("Networks Header Files" FILES ${NN_INCLUDE_NETWORKS})
source_group("Optimizers Header Files" FILES ${NN_INCLUDE_OPTIMIZERS})
source_group("Source Files" FILES ${NN_SRC})

add_library(NN
	${NN_INCLUDE_NETWORKS}
	${NN_INCLUDE_OPTIMIZERS}
	${NN_SRC}
)

target_include_directories(NN
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/Networks"
	PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/Optimizers"
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(NN PRIVATE Eigen3::Eigen)