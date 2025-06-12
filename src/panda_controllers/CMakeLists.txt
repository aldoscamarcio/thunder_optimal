cmake_minimum_required(VERSION 2.8.3)
project(panda_controllers)

set(CMAKE_BUILD_TYPE release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(catkin REQUIRED COMPONENTS
	controller_interface
	dynamic_reconfigure
	franka_hw
	geometry_msgs
	franka_msgs
	hardware_interface
	# message_generation #
	pluginlib
	realtime_tools
	roscpp
	rospy
	sensor_msgs
	visualization_msgs
)

find_package(Eigen3 REQUIRED)
find_package(Franka 0.9.2 REQUIRED)
find_package(NLopt REQUIRED)
find_package(yaml-cpp REQUIRED)

# Stampa il percorso alle directory di inclusione della libreria Nlopt
# message(STATUS "NLopt include directories: ${EIGEN3_INCLUDE_DIRS}")

# Stampa le librerie Nlopt
# message(STATUS "NLopt libraries: ${NLOPT_LIBRARIES}")
# find_package(nlopt REQUIRED)
#find_package(Boost REQUIRED)

#find_package(CasADi REQUIRED) # for regressor
#find_library(CASADI_LIB casadi)

# # to add custom msg
# add_message_files(FILES
# 	desTrajEE.msg
# 	log_adaptive_cartesian.msg
# 	log_adaptive_joints.msg
# 	point.msg
# 	link_params.msg
# 	flag.msg
# 	Vec7D.msg
# 	rpy.msg
# 	impedanceGain.msg
# )
## to create msg, specify dependecies
# generate_messages(
#   DEPENDENCIES
#   std_msgs
#   geometry_msgs
# )

include_directories(
	include #
	${catkin_INCLUDE_DIRS} 
	${EIGEN3_INCLUDE_DIRS}
	${NLOPT_INCLUDE_DIRS}
	${YAML_CPP_INCLUDE_DIRS}  # Include yaml-cpp directories
	# ${Boost_INCLUDE_DIRS}
	# ${CasADi_INCLUDE_DIRS}
	# ${CMAKE_CURRENT_SOURCE_DIR}/include #
)


catkin_package(
	INCLUDE_DIRS include #
	LIBRARIES ${PROJECT_NAME}
	CATKIN_DEPENDS
	controller_interface
	dynamic_reconfigure
	franka_hw
	geometry_msgs
	franka_msgs
	hardware_interface
	message_runtime
	pluginlib
	realtime_tools
	roscpp
	# DEPENDS Franka #
)

add_library(${PROJECT_NAME}   
	src/computed_torque.cpp 	
	src/franka_gen.cpp
	src/thunder_franka.cpp
	# src/optimal_control.cpp
	src/thunder_optimization.cpp
)

target_link_libraries(${PROJECT_NAME} PUBLIC
	# ${Franka_LIBRARIES} #
	${catkin_LIBRARIES}
	yaml-cpp
	nlopt
	#${Boost_LIBRARIES}
	#${CASADI_LIB}
)

# add_dependencies(${PROJECT_NAME}
# 	${${PROJECT_NAME}_EXPORTED_TARGETS}
# 	panda_controllers ${catkin_EXPORTED_TARGETS}
# )


add_executable(menu
 	src/menu.cpp
	# src/thunder_franka.cpp
	# src/franka_gen.cpp
)

target_link_libraries(menu
	${PROJECT_NAME}
	${catkin_LIBRARIES}
	${NLOPT_LIBRARIES}
	yaml-cpp
)

# add_dependencies(menu
# 	${${PROJECT_NAME}_EXPORTED_TARGETS} 
# 	${catkin_EXPORTED_TARGETS}
# )

# target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
# 	${Franka_INCLUDE_DIRS}
# 	${EIGEN3_INCLUDE_DIRS}
# 	${catkin_INCLUDE_DIRS}
# 	# ${nlopt_INCLUDE_DIRS}
# )
# target_include_directories(${PROJECT_NAME} PUBLIC
# 	include
# )

## Installation
# install(TARGETS ${PROJECT_NAME}
# 	ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
# 	LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
# 	RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )
# install(DIRECTORY launch
# 	DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )
# install(DIRECTORY config
# 	DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )
# install(FILES ${PROJECT_NAME}_plugin.xml
# 	DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )