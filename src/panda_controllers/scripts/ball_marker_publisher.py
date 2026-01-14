#!/usr/bin/env python  
import rospy
from visualization_msgs.msg import Marker

def ball_publisher():   
    rospy.init_node('ball_marker_publisher')
    pub = rospy.Publisher('/ball_marker', Marker, queue_size=10)

    marker = Marker()
    marker.header.frame_id = "panda_link0"
    marker.header.stamp = rospy.Time.now()
    marker.ns = "obstacle"
    marker.id = 1
    marker.type = Marker.SPHERE
    marker.action = Marker.ADD
    marker.pose.position.x = -0.15126  
    marker.pose.position.y = 0.336159
    marker.pose.position.z = 0.545768
    marker.pose.orientation.w = 1.0
    marker.scale.x = 0.1  # diametro della sfera
    marker.scale.y = 0.1
    marker.scale.z = 0.1
    marker.color.r = 0.0
    marker.color.g = 1.0
    marker.color.b = 0.0
    marker.color.a = 0.8

    rate = rospy.Rate(1)  # Pubblica ogni secondo
    while not rospy.is_shutdown():
        marker.header.stamp = rospy.Time.now()
        pub.publish(marker)
        rate.sleep()

if __name__ == '__main__':
    ball_publisher()
