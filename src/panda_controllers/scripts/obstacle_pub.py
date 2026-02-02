#!/usr/bin/env python3
import rospy
from visualization_msgs.msg import Marker, MarkerArray

# Funzione helper aggiornata per accettare frame_id
def create_marker(id, marker_type, pose, size, color, frame_id):
    marker = Marker()
    marker.header.frame_id = frame_id
    marker.header.stamp = rospy.Time.now()
    marker.ns = "gazebo_obstacles"
    marker.id = id
    marker.type = marker_type
    marker.action = Marker.ADD
    
    # Posizione
    marker.pose.position.x = pose[0]
    marker.pose.position.y = pose[1]
    marker.pose.position.z = pose[2]
    marker.pose.orientation.w = 1.0 
    
    # Dimensione
    marker.scale.x = size[0]
    marker.scale.y = size[1]
    marker.scale.z = size[2]
    
    # Colore
    marker.color.r = color[0]
    marker.color.g = color[1]
    marker.color.b = color[2]
    marker.color.a = color[3]
    
    return marker

def obstacle_publisher():
    rospy.init_node('obstacle_visualizer_node')
    
    # Legge il parametro frame_id (o usa 'world' se non lo trova)
    frame_id_param = rospy.get_param("~frame_id", "world")
    
    pub = rospy.Publisher('/Obstacle_Marker_Array', MarkerArray, queue_size=10)
    rate = rospy.Rate(1) # 1 Hz è sufficiente per ostacoli statici

    while not rospy.is_shutdown():
        marker_array = MarkerArray()

        # 1. Sfera 1
        # s1 = create_marker(1, Marker.SPHERE, 
        #                    [0.40, 0, 0.65], 
        #                    [0.14, 0.14, 0.14], 
        #                    [0.0, 0.0, 1.0, 1.0],
        #                    frame_id_param)
        
        # 2. Muro
        w1 = create_marker(2, Marker.CUBE, 
                           [0.35, 0, 0.60], 
                           [0.2, 0.02, 0.2], 
                           [1.0, 0.0, 0.0, 1.0],
                           frame_id_param)
        
        # # 3. Sfera 2
        # s2 = create_marker(3, Marker.SPHERE, 
        #                    [0.30, 0.45, 0.80], 
        #                    [0.16, 0.16, 0.16], 
        #                    [0.0, 0.0, 1.0, 1.0],
        #                    frame_id_param)

        # marker_array.markers.append(s1)
        marker_array.markers.append(w1)
        # marker_array.markers.append(s2)

        pub.publish(marker_array)
        rate.sleep()

if __name__ == '__main__':
    try:
        obstacle_publisher()
    except rospy.ROSInterruptException:
        pass