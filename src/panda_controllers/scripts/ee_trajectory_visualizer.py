#!/usr/bin/env python3
import rospy
import tf
import traceback 
from visualization_msgs.msg import Marker
from geometry_msgs.msg import Point
from std_msgs.msg import ColorRGBA

class TrajectoryVisualizer:
    def __init__(self):
        rospy.init_node('ee_trajectory_visualizer')
        
        self.listener = tf.TransformListener()
        
        self.line_pub = rospy.Publisher('/ee_vis/line', Marker, queue_size=10)
        self.points_pub = rospy.Publisher('/ee_vis/points', Marker, queue_size=10)
        
        self.trajectory_points = []
        
        # Memoria di 5000 punti (circa 4 minuti)
        self.max_points = 5000 
        
        rospy.sleep(1.0)
        self.rate = rospy.Rate(20.0) 
        self.run()

    def run(self):
        rospy.loginfo("Trajectory Visualizer Started.")
        
        while not rospy.is_shutdown():
            try:
                (trans, rot) = self.listener.lookupTransform('panda_link0', 'panda_EE', rospy.Time(0))
                self.add_point(trans)
                self.publish_visuals()
            except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
                pass
            except Exception as e:
                rospy.logerr("Errore: " + str(e))
                
            self.rate.sleep()

    def add_point(self, trans):
        p = Point()
        p.x = trans[0]
        p.y = trans[1]
        p.z = trans[2]
        
        # Ottimizzazione: se il robot è fermo (spostamento < 2mm), non aggiungiamo punti
        if self.trajectory_points:
            last_p = self.trajectory_points[-1]
            dist = ((p.x - last_p.x)**2 + (p.y - last_p.y)**2 + (p.z - last_p.z)**2)**0.5
            if dist < 0.002: 
                return

        self.trajectory_points.append(p)
        
        if len(self.trajectory_points) > self.max_points:
            self.trajectory_points.pop(0)

    def publish_visuals(self):
        if not self.trajectory_points:
            return

        current_time = rospy.Time.now()

        # --- MARKER 1: LINEA SOLIDA ---
        line_marker = Marker()
        line_marker.header.frame_id = "panda_link0"
        line_marker.header.stamp = current_time
        line_marker.ns = "ee_trail"
        line_marker.id = 0
        line_marker.type = Marker.LINE_STRIP
        line_marker.action = Marker.ADD
        line_marker.pose.orientation.w = 1.0
        
        # SPESSORE DELLA LINEA (Aumentato)
        line_marker.scale.x = 0.01  # 1 cm di larghezza (era 0.005)
        
        # COLORE UNICO (Verde)
        line_marker.color.r = 0.0
        line_marker.color.g = 1.0
        line_marker.color.b = 0.0
        line_marker.color.a = 1.0  # Opacità completa (no trasparenza)
        
        line_marker.points = self.trajectory_points
        self.line_pub.publish(line_marker)

        # --- MARKER 2: PUNTI (SFERE) ---
        points_marker = Marker()
        points_marker.header.frame_id = "panda_link0"
        points_marker.header.stamp = current_time
        points_marker.ns = "ee_dots"
        points_marker.id = 1
        points_marker.type = Marker.SPHERE_LIST
        points_marker.action = Marker.ADD
        points_marker.pose.orientation.w = 1.0
        
        # DIMENSIONE PUNTI (Leggermente più grandi della linea)
        points_marker.scale.x = 0.02
        points_marker.scale.y = 0.02
        points_marker.scale.z = 0.02
        
        # Colore uguale alla linea
        points_marker.color = line_marker.color
        
        # Mostriamo 1 punto ogni 10 per leggerezza, ma sembreranno connessi dalla linea
        points_marker.points = self.trajectory_points[::10] 
        
        self.points_pub.publish(points_marker)

if __name__ == '__main__':
    try:
        TrajectoryVisualizer()
    except rospy.ROSInterruptException:
        pass