# generated from catkin/cmake/template/pkg.context.pc.in
CATKIN_PACKAGE_PREFIX = ""
PROJECT_PKG_CONFIG_INCLUDE_DIRS = "${prefix}/include".split(';') if "${prefix}/include" != "" else []
PROJECT_CATKIN_DEPENDS = "controller_interface;dynamic_reconfigure;franka_hw;geometry_msgs;franka_msgs;hardware_interface;message_runtime;pluginlib;realtime_tools;roscpp".replace(';', ' ')
PKG_CONFIG_LIBRARIES_WITH_PREFIX = "-lpanda_controllers".split(';') if "-lpanda_controllers" != "" else []
PROJECT_NAME = "panda_controllers"
PROJECT_SPACE_DIR = "/home/franko/Scrivania/thunder_optimal/install"
PROJECT_VERSION = "1.0.0"
