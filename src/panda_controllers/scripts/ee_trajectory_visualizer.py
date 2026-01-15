#!/usr/bin/env python  
import rospy
import tf
from geometry_msgs.msg import Point, PoseStamped, Quaternion
from visualization_msgs.msg import Marker
from nav_msgs.msg import Path

class TrajectoryVisualizer:
    def __init__(self):
        rospy.init_node('ee_trajectory_visualizer')
        self.listener = tf.TransformListener()
        self.marker_pub = rospy.Publisher('/ee_marker_trajectory', Marker, queue_size=10)
        self.path_pub = rospy.Publisher('/ee_path_trajectory', Path, queue_size=10)
        self.trajectory_points = []
        self.path_msg = Path()
        self.path_msg.header.frame_id = "panda_link0"
        self.max_points = 5000
        self.rate = rospy.Rate(20.0)
        self.run()

    def run(self):
        while not rospy.is_shutdown():
            try:
                (trans, rot) = self.listener.lookupTransform('panda_link0', 'panda_EE', rospy.Time(0))
                self.add_point(trans, rot)
                self.publish_marker()
                self.publish_path()
            except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
                pass
            self.rate.sleep()

    def add_point(self, trans, rot):
        p = Point(trans[0], trans[1], trans[2])
        self.trajectory_points.append(p)
        if len(self.trajectory_points) > self.max_points:
            self.trajectory_points.pop(0)

        pose = PoseStamped()
        pose.header.stamp = rospy.Time.now()
        pose.header.frame_id = "panda_link0"
        pose.pose.position.x, pose.pose.position.y, pose.pose.position.z = trans
        pose.pose.orientation.x, pose.pose.orientation.y, pose.pose.orientation.z, pose.pose.orientation.w = rot
        self.path_msg.poses.append(pose)
        if len(self.path_msg.poses) > self.max_points:
            self.path_msg.poses.pop(0)

    def publish_marker(self):
        marker = Marker()
        marker.header.frame_id = "panda_link0"
        marker.header.stamp = rospy.Time.now()
        marker.ns = "ee_trajectory"
        marker.id = 0
        marker.type = Marker.LINE_STRIP
        marker.action = Marker.ADD
        
        marker.pose.position.x = 0.0
        marker.pose.position.y = 0.0
        marker.pose.position.z = 0.0
        marker.pose.orientation.x = 0.0
        marker.pose.orientation.y = 0.0
        marker.pose.orientation.z = 0.0
        marker.pose.orientation.w = 1.0  # Quaternione identità
        
        # PROPRIETÀ VISIVE
        marker.scale.x = 0.02  # Linea spessa
        marker.color.r = 1.0   # Rosso puro
        marker.color.g = 0.0
        marker.color.b = 0.0
        marker.color.a = 1.0   # Opaco
        
        # Altri parametri opzionali
        marker.lifetime = rospy.Duration()  # Marker permanente
        marker.frame_locked = True  # Bloccato al frame (si muove con il robot)
        
        marker.points = self.trajectory_points
        self.marker_pub.publish(marker)

    def publish_path(self):
        self.path_msg.header.stamp = rospy.Time.now()
        self.path_pub.publish(self.path_msg)

if __name__ == '__main__':
    TrajectoryVisualizer()