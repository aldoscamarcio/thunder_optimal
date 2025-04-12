#include <iostream>
#include <eigen3/Eigen/Dense>

#include <unistd.h>
#include <cstdlib>
#include <signal.h>
#include <math.h>
#include <vector>
#include <geometry_msgs/PoseStamped.h>
#include "ros/ros.h"
#include "utils/thunder_franka.h"
#include "utils/thunder_optimization.h"

#include <sstream>
// #include <eigen_conversions/eigen_msg.h>
#include "nlopt.hpp"
// ROS Service and Message Includes
#include "std_msgs/Float64.h"
#include "std_msgs/Bool.h"
#include "std_srvs/SetBool.h"
// #include "geometry_msgs/Pose.h"
#include "sensor_msgs/JointState.h"

const std::string conf_file = "../config/franka_conf.yaml";

using namespace std;
using std::cout;
using std::endl;

#define alpha 0.1
bool init_flag = false;
bool init_q0 = false;

struct traj_struct
{
	Eigen::Matrix<double, 7, 1> pos_des;
	Eigen::Matrix<double, 7, 1> vel_des;
	Eigen::Matrix<double, 7, 1> acc_des;
} traj;

// define q0 as 7x1 matrix
Eigen::Matrix<double, 7, 1> q0;
const double q_lim_upp[] = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
const double q_lim_low[] = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};

// Define the function toDouble() be called when ctrl-c (SIGINT) is sent toDouble() process
void signal_callback_handler(int signum)
{
	cout << "Caught signal " << signum << endl;
	// Terminate program
	exit(signum);
}

// void poseCallback(const geometry_msgs::PoseStampedConstPtr& msg) {
//   pos << msg->pose.position.x, msg->pose.position.y, msg->pose.position.z;
//   orient << msg->pose.orientation.x, msg->pose.orientation.y, msg->pose.orientation.z;
// }

void jointsCallback(const sensor_msgs::JointStateConstPtr &msg)
{
	q0 = Eigen::Map<const Eigen::Matrix<double, 7, 1>>((msg->position).data());
	init_q0 = true;
}

void interpolator_pos(Eigen::Matrix<double, 7, 1> pos_i, Eigen::Matrix<double, 7, 1> pos_f, double tf, double t)
{
	traj.pos_des << pos_i + (pos_i - pos_f) * (15 * pow((t / tf), 4) - 6 * pow((t / tf), 5) - 10 * pow((t / tf), 3));
	traj.vel_des << (pos_i - pos_f) * (60 * (pow(t, 3) / pow(tf, 4)) - 30 * (pow(t, 4) / pow(tf, 5)) - 30 * (pow(t, 2) / pow(tf, 3)));
	traj.acc_des << (pos_i - pos_f) * (180 * (pow(t, 2) / pow(tf, 4)) - 120 * (pow(t, 3) / pow(tf, 5)) - 60 * (t / pow(tf, 3)));
}

// void demo_inf_XY(Eigen::Vector3d pos_i, double t){
// 	Eigen::Vector3d tmp;
// 	tmp << sin(t)/8, sin(t/2)/4, 0;
// 	traj.pos_des << pos_i + tmp;
// 	traj.vel_des << cos(t)/8, cos(t/2)/8, 0;
// 	traj.acc_des << -sin(t)/8, -sin(t/2)/16, 0;
// }

// void demo_inf_XYZ(Eigen::Vector3d pos_i, double t,double zf,double tf){
// 	Eigen::Vector3d tmp;
// 	tmp << sin(t)/8, sin(t/2)/4, ((zf-pos_i(2))/tf)*t;
// 	traj.pos_des << pos_i + tmp;
// 	traj.vel_des << cos(t)/8, cos(t/2)/8, (zf-pos_i(2))/tf;
// 	traj.acc_des << -sin(t)/8, -sin(t/2)/16, 0;
// }

// void demo_circle_xy(Eigen::Vector3d pos_i, double t,double zf,double tf){
//   Eigen::Vector3d tmp;
//   tmp << 0.1*cos(t), 0.1*sin(t), ((zf-pos_i(2))/tf)*t;
//   traj.pos_des << pos_i + tmp;
//   traj.vel_des << -0.1*sin(t), 0.1*cos(t), (zf-pos_i(2))/tf;
//   traj.acc_des << -0.1*cos(t), -0.1*sin(t), 0;
// }

int main(int argc, char **argv)
{
	ros::init(argc, argv, "menu");

	ros::NodeHandle node_handle;
	thunder_franka robot;
	robot.load_conf(conf_file);
	int NJ = robot.get_numJoints();
	cout << "NJ" << NJ << endl;
	Eigen::VectorXd q(NJ), dq(NJ), dqr(NJ), ddqr(NJ);

	/* Test */
	q = q.setOnes();
	dq = dq.setOnes();
	dqr = dqr.setOnes();
	ddqr = ddqr.setOnes();

	robot.setArguments(q, dq, dqr, ddqr);

	ros::Publisher pub_cmd = node_handle.advertise<sensor_msgs::JointState>("/computed_torque_controller/command", 1000);
	ros::Subscriber sub_joints = node_handle.subscribe<sensor_msgs::JointState>("/franka_state_controller/joint_states", 1, &jointsCallback);
	// ros::Subscriber sub_pose =  node_handle.subscribe("/franka_state_controller/franka_ee_pose", 1, &poseCallback);

	// creating trajectory message
	sensor_msgs::JointState traj_msg;

	// SET SLEEP TIME 1000 ---> 1 kHz
	ros::Rate loop_rate(10);

	srand(time(NULL));
	double tf;
	Eigen::Matrix<double, 7, 1> q_int;
	Eigen::Matrix<double, 7, 1> qf;
	Eigen::Matrix<double, 7, 1> v0;
	Eigen::Matrix<double, 7, 1> vf;
	Eigen::Matrix<double, 7, 1> a0;
	Eigen::Matrix<double, 7, 1> af;

	v0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	a0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	int campioni = 101;
	//Eigen::VectorXd POS_INIT(NJ * campioni), VEL_INIT(NJ * campioni), ACC_INIT(NJ * campioni);
	Eigen::VectorXd POS_INIT(NJ), VEL_INIT(NJ), ACC_INIT(NJ);
	int time_step = 0;
	XmlRpc::XmlRpcValue menu_par;

	// Initialize Ctrl-C
	signal(SIGINT, signal_callback_handler);

	ros::Time t_init;
	init_q0 = false;
	double t = 0;
	int choice;
	int demo = -1;
	int yaml = 0;

	while (ros::ok())
	{
		demo = -1;
		if (yaml == 1)
		{
			choice = 5;
		}
		else
		{
			cout << "choice:   (1: joints min-jerk,  2: go toDouble() init,  3: go toDouble() random  4: yaml 6: optimal_min-jerk) " << endl;
			cin >> choice;
		}
		if (choice == 1)
		{
			cout << "duration: " << endl;
			cin >> tf;
			cout << "final_joint_positions: " << endl;
			cin >> qf(0);
			cin >> qf(1);
			cin >> qf(2);
			cin >> qf(3);
			cin >> qf(4);
			cin >> qf(5);
			cin >> qf(6);
		}
		else if (choice == 2)
		{
			std::vector<double> qf_array;
			if (!node_handle.getParam("/menu/Q0_INIT", qf_array))
				ROS_ERROR("Failed toDouble() get parameter from server.");
			qf = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(qf_array.data(), qf_array.size());
			choice = 1;
			tf = 3.0;
		}
		else if (choice == 3)
		{
			for (int i = 0; i < 7; i++)
			{
				double q_low = q_lim_low[i];
				double q_upp = q_lim_upp[i];
				qf(i) = q_low + (float(rand()) / RAND_MAX) * (q_upp - q_low);
			}
			choice = 1;
			tf = 3.0;
		}
		else if (choice == 4)
		{
			cout << "-not implemented yet-" << endl;
		}
		else if (choice == 6)
		{
			tf = 10;

			q_int << -1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972; // Posizioni iniziali

			qf << 0.0, 0.0, 0.0, -0.2, 2.0, 1.0, 0.0; // Posizioni finali

			v0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità iniziali

			vf << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità finali

			a0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni iniziali

			af << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni finali

			// cout << "duration: " << endl;
			// cin >> tf;
			// cout << "final_joint_positions: " << endl;
			// cin >> qf(0);
			// cin >> qf(1);
			// cin >> qf(2);
			// cin >> qf(3);
			// cin >> qf(4);
			// cin >> qf(5);
			// cin >> qf(6);
			// cin >> vf(0);
			// cin >> vf(1);
			// cin >> vf(2);
			// cin >> vf(3);
			// cin >> vf(4);
			// cin >> vf(5);
			// cin >> vf(6);
			// cin >> af(0);
			// cin >> af(1);
			// cin >> af(2);
			// cin >> af(3);
			// cin >> af(4);
			// cin >> af(5);
			// cin >> af(6);
		}
		if (!init_flag)
		{
			// cout << "frankino" << endl;
			// cout << "q0"<< q_int<<endl;
			init_q0 = false;
			init_flag = true;
		}

		ros::spinOnce();

		t_init = ros::Time::now();

		t = (ros::Time::now() - t_init).toSec();

		while (t <= tf && init_q0)
		{
			if (choice == 1)
			{
				interpolator_pos(q0, qf, tf, t);
			}
			else if (choice == 4)
			{
				break;
			}
			else
			{
				break;
			}

			traj_msg.header.stamp = ros::Time::now();

			std::vector<double> pos_des{traj.pos_des[0], traj.pos_des[1], traj.pos_des[2], traj.pos_des[3], traj.pos_des[4], traj.pos_des[5], traj.pos_des[6]};
			traj_msg.position = pos_des;
			std::vector<double> vel_des{traj.vel_des[0], traj.vel_des[1], traj.vel_des[2], traj.vel_des[3], traj.vel_des[4], traj.vel_des[5], traj.vel_des[6]};
			traj_msg.velocity = vel_des;
			std::vector<double> acc_des{traj.acc_des[0], traj.acc_des[1], traj.acc_des[2], traj.acc_des[3], traj.acc_des[4], traj.acc_des[5], traj.acc_des[6]};
			traj_msg.effort = acc_des;
			pub_cmd.publish(traj_msg);

			loop_rate.sleep();

			t = (ros::Time::now() - t_init).toSec();
		}

		if (choice == 6)
		{
			// Calcolo dei coefficienti per ciascun giunto
			std::vector<std::vector<double>> joint_coeffs(NJ);
			for (int i = 0; i < NJ; ++i)
			{
				joint_coeffs[i] = calculateCoefficients(q_int[i], qf[i], v0[i], vf[i], a0[i], af[i], t_init.toSec(), tf);
				std::cout << "coeffs: " << joint_coeffs[i][0] << ", " << joint_coeffs[i][1] << ", " << joint_coeffs[i][2] << ", " << joint_coeffs[i][3] << ", " << joint_coeffs[i][4] << ", " << joint_coeffs[i][5] << std::endl;


			}
			while (t <= tf)
			{
				std::cout << "Time: " << t << std::endl;
				// Generazione della traiettoria
				for (int i = 0; i < NJ; i++)
				{
					double pos, vel, acc, tau;
					calculateTrajectory(t, t_init.toSec(), joint_coeffs[i], pos, vel, acc);
					std::cout << "  Joint " << i << ": Pos = " << pos << ", Vel = " << vel << ", Acc = " << acc << std::endl;
					// POS(i, time_step) = pos;
					// VEL(i, time_step) = vel;
					// ACC(i, time_step) = acc;

					// // Incolonna tutte le varibili di giunto ad ogni step
					// POS_INIT(time_step * NJ + i) = pos;
					// VEL_INIT(time_step * NJ + i) = vel;
					// ACC_INIT(time_step * NJ + i) = acc;
					POS_INIT(i)=pos;
					VEL_INIT(i)=vel;
					ACC_INIT(i)=acc;
				}

				time_step++;
				std::vector<double> pos_des{POS_INIT[0],POS_INIT[1],POS_INIT[2],POS_INIT[3],POS_INIT[4],POS_INIT[5],POS_INIT[6]};
				traj_msg.position = pos_des;
				std::vector<double> vel_des{VEL_INIT[0],VEL_INIT[1],VEL_INIT[2],VEL_INIT[3],VEL_INIT[4],VEL_INIT[5],VEL_INIT[6]};
				traj_msg.velocity = vel_des;
				std::vector<double> acc_des{ACC_INIT[0],ACC_INIT[1],ACC_INIT[2],ACC_INIT[3],ACC_INIT[4],ACC_INIT[5],ACC_INIT[6]};
				traj_msg.effort = acc_des;
				pub_cmd.publish(traj_msg);

				loop_rate.sleep();

				t = (ros::Time::now() - t_init).toSec();

				// cout<<POS_INIT<<endl;
			}

			std::cout << "time_step:" << time_step << std::endl;
			// for (int i = 0; i < NJ * campioni; i++)
			// {
			// 	std::cout << "POS_INIT: " << i << ": " << POS_INIT[i] << endl;
			// }
		}
	}
	return 0;
}