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
	//cout << "NJ" << NJ << endl;
	Eigen::VectorXd q(NJ), dq(NJ), dqr(NJ), ddqr(NJ);

	ros::Publisher pub_cmd = node_handle.advertise<sensor_msgs::JointState>("/computed_torque_controller/command", 1000);
	ros::Subscriber sub_joints = node_handle.subscribe<sensor_msgs::JointState>("/franka_state_controller/joint_states", 1, &jointsCallback);
	// ros::Subscriber sub_pose =  node_handle.subscribe("/franka_state_controller/franka_ee_pose", 1, &poseCallback);

	// creating trajectory message
	sensor_msgs::JointState traj_msg;

	// SET SLEEP TIME 1000 ---> 1 kHz
	int frequenza=10; // Hz
	ros::Rate loop_rate(frequenza);

	srand(time(NULL));
	double tf;
	Eigen::Matrix<double, 7, 1> q_int;
	Eigen::Matrix<double, 7, 1> qf;
	Eigen::Matrix<double, 7, 1> v0;
	Eigen::Matrix<double, 7, 1> vf;
	Eigen::Matrix<double, 7, 1> a0;
	Eigen::Matrix<double, 7, 1> af;

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
		int time_step = 0;
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

			//q_int << -1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972; // Posizioni iniziali

			qf << -1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972; // Posizioni finali

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

		int campioni = tf*frequenza+1; // Numero di campioni
		int size_q = NJ * campioni; // Dimensione di q
		Eigen::VectorXd POS_INIT(NJ * campioni), VEL_INIT(NJ * campioni), ACC_INIT(NJ * campioni);

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

		if (choice == 6 && init_q0)
		{
			q_int=q0;

			// cout << "q0 "<< q0 << endl;
			// cout << "q_int " << q_int << endl;

			// Calcolo dei coefficienti per ciascun giunto
			std::vector<std::vector<double>> joint_coeffs(NJ);
			for (int i = 0; i < NJ; ++i)
			{
				joint_coeffs[i] = calculateCoefficients(q_int[i], qf[i], v0[i], vf[i], a0[i], af[i], t_init.toSec(), tf);
				//std::cout << "coeffs: " << joint_coeffs[i][0] << ", " << joint_coeffs[i][1] << ", " << joint_coeffs[i][2] << ", " << joint_coeffs[i][3] << ", " << joint_coeffs[i][4] << ", " << joint_coeffs[i][5] << std::endl;
			}
			while (t <= tf)
			{
				// std::cout << "Time: " << t << std::endl;
				// Generazione della traiettoria
				for (int i = 0; i < NJ; i++)
				{

					double pos, vel, acc;
					calculateTrajectory(t, t_init.toSec(), joint_coeffs[i], pos, vel, acc);
					//std::cout << "  Joint " << i << ": Pos = " << pos << ", Vel = " << vel << ", Acc = " << acc << std::endl;
					// // In riga tutti i giunti e colonn i vari step
					// POS(i, time_step) = pos;
					// VEL(i, time_step) = vel;
					// ACC(i, time_step) = acc;

					// Incolonna e colleziona tutte le varibili di giunto per ogni step
					POS_INIT(time_step * NJ + i) = pos;
					VEL_INIT(time_step * NJ + i) = vel;
					ACC_INIT(time_step * NJ + i) = acc;

					// // Incolonna tutte le varibili di giunto ad ogni step
					// POS_INIT(i) = pos;
					// VEL_INIT(i) = vel;
					// ACC_INIT(i) = acc;
				}

				

				// std::vector<double> pos_des{POS_INIT[time_step * NJ + 0], POS_INIT[time_step * NJ +1], POS_INIT[time_step * NJ +2], POS_INIT[time_step * NJ +3], POS_INIT[time_step * NJ +4], POS_INIT[time_step * NJ +5], POS_INIT[time_step * NJ +6]};
				// traj_msg.position = pos_des;
				// std::vector<double> vel_des{VEL_INIT[time_step * NJ +0], VEL_INIT[time_step * NJ +1], VEL_INIT[time_step * NJ +2], VEL_INIT[time_step * NJ +3], VEL_INIT[time_step * NJ +4], VEL_INIT[time_step * NJ +5], VEL_INIT[time_step * NJ +6]};
				// traj_msg.velocity = vel_des;
				// std::vector<double> acc_des{ACC_INIT[time_step * NJ +0], ACC_INIT[time_step * NJ +1], ACC_INIT[time_step * NJ +2], ACC_INIT[time_step * NJ +3], ACC_INIT[time_step * NJ +4], ACC_INIT[time_step * NJ +5], ACC_INIT[time_step * NJ +6]};
				// traj_msg.effort = acc_des;
				// pub_cmd.publish(traj_msg);
				
				time_step++;

				loop_rate.sleep();

				t = (ros::Time::now() - t_init).toSec();

				
			}


			OptimizationData optData;
			optData.robot; // Inizializza l'oggetto robot
			optData.q0 = q_int;
			optData.v0 = v0;
			optData.a0 = a0;
			optData.qf = qf;
			optData.vf = vf;
			optData.af = af;
			optData.size_q = size_q;
			optData.dt = 1.0/frequenza; // Passo temporale
			optData.campioni = campioni;

			// define the optimization problem
			nlopt::opt opt(nlopt::LD_MMA, 3 * NJ * campioni);
			opt.set_min_objective(objective, &optData);
			opt.set_xtol_rel(1e-3);

			std::vector<double> ub(3 * NJ * campioni), lb(3 * NJ * campioni), ubq(NJ), lbq(NJ), ubdq(NJ), lbdq(NJ), ubddq(NJ), lbddq(NJ);
			lbq = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
			ubq = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
			lbdq = {-2.175, -2.175, -2.175, -2.175, -2.61, -2.61, -2.61};
			ubdq = {2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61};
			lbddq = {-15, -7.5, -10, -12.5, -15, -20, -20};
			ubddq = {15, 7.5, 10, 12.5, 15, 20, 20};

			float tol = 4e-2;
			int j=0;
			// Vincoli sulle condizioni iniziali
			if (j == 0)
			{
				for (int i = 0; i < NJ; i++)
				{
					lb[j * NJ + i] = q_int[i] - tol;
					ub[j * NJ + i] = q_int[i] + tol;
					lb[j * NJ + size_q + i] = v0[i] - tol;
					ub[j * NJ + size_q + i] = v0[i] + tol;
					lb[j * NJ + 2 * size_q + i] = a0[i] - tol;
					ub[j * NJ + 2 * size_q + i] = a0[i] + tol;
				}
			}

			// Vincoli sulle condizioni intermedie

			for (int j = 1; j < campioni-1; j++)
			{
				for (int i = 0; i < NJ; i++)
				{
					lb[j * NJ + i] = lbq[i];
					ub[j * NJ + i] = ubq[i];
					lb[j * NJ + size_q + i] = lbdq[i];
					ub[j * NJ + size_q + i] = ubdq[i];
					lb[j * NJ + 2 * size_q + i] = lbddq[i];
					ub[j * NJ + 2 * size_q + i] = ubddq[i];
				}
			}

			j=campioni-1;

			// Vincoli sulle condizioni finali

			if(j==campioni-1)
			{
				for (int i = 0; i < NJ; i++)
				{
					lb[j * NJ + i] = qf[i] - tol;
					//std::cout << "lb " << j * NJ + i << ":" << lb[j * NJ + i] << std::endl;
					ub[j * NJ + i] = qf[i] + tol;
					//std::cout << "ub " << j * NJ + i << ":" << ub[j * NJ + i] << std::endl;
					lb[j * NJ + size_q + i] = vf[i] - tol;
					//std::cout << "lb " << j * NJ + size_q + i << ":" << lb[j * NJ + size_q + i] << std::endl;
					ub[j * NJ + size_q + i] = vf[i] + tol;
					//std::cout << "ub " << j * NJ + size_q + i << ":" << ub[j * NJ + size_q + i] << std::endl;
					lb[j * NJ + 2 * size_q + i] = af[i] - tol;
					//std::cout << "lb " << j * NJ + 2 * size_q + i << ":" << lb[j * NJ + 2 * size_q + i] << std::endl;
					ub[j * NJ + 2 * size_q + i] = af[i] + tol;
					//std::cout << "ub " << j * NJ + 2 * size_q + i << ":" << ub[j * NJ + 2 * size_q + i] << std::endl;
				}
			}

			// for (int i = 0; i < lb.size(); i++)
			// {
			// 	std::cout << "lb" << i << ": " << lb[i]<<" ---- " << "ub" << i <<": " << ub[i] << std::endl;
			// }
			
			opt.set_upper_bounds(ub);
			opt.set_lower_bounds(lb);
			
			std::vector<ConsistencyConstraintIneq> constraints;
			double dt=1.0/frequenza;
			const double eps = 4e-2;

			int numero_totale_vincoli = (campioni - 1) * 6; // <-- Calcola il numero totale

			if (numero_totale_vincoli > 0) { //evita di chiamare reserve(0)
				constraints.reserve(numero_totale_vincoli);
			}

			for (int k = 0; k < campioni-1; k++) {
				// Posizione
				constraints.push_back({k, NJ, size_q, dt, 0, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			
				constraints.push_back({k, NJ, size_q, dt, 0, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			
				// Velocità
				constraints.push_back({k, NJ, size_q, dt, 1, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			
				constraints.push_back({k, NJ, size_q, dt, 1, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			
				// Accelerazione
				constraints.push_back({k, NJ, size_q, dt, 2, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			
				constraints.push_back({k, NJ, size_q, dt, 2, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);
			}
			

			// define the initial guess
			std::vector<double> vettore(3 * NJ * campioni); // Inizializza il vettore x
			for (int i = 0; i < POS_INIT.size(); i++)
			{
				vettore[i] = POS_INIT[i];
				vettore[i + POS_INIT.size()] = VEL_INIT[i];
				vettore[i + POS_INIT.size() + VEL_INIT.size()] = ACC_INIT[i];
			}

			double minf;
			nlopt::result result = opt.optimize(vettore, minf);

			// Print the optimized values
			printf("Optimized values:\n");
			printf("q_opt:\n");
			for (int i = 0; i < 3*POS_INIT.size(); i++)
			{	
				if(i== POS_INIT.size())
				{
					printf("\n dq_opt:\n");
				}
				if(i== 2*POS_INIT.size())
				{
					printf("\nddq_opt:\n");
				}
				printf("%g ", vettore[i]);
				printf(",\n");

			}


			std::cout << "time_step:" << time_step << std::endl;
			// for (int i = 0; i < NJ * campioni; i++)
			// {
			// 	std::cout << "POS_INIT: " << i << ": " << POS_INIT[i] << endl;
			// }

			time_step = 0;
		while (time_step < campioni)
		{
			std::vector<double> pos_des{vettore[time_step * NJ + 0], vettore[time_step * NJ +1], vettore[time_step * NJ +2], vettore[time_step * NJ +3], vettore[time_step * NJ +4], vettore[time_step * NJ +5], vettore[time_step * NJ +6]};
				traj_msg.position = pos_des;
				std::vector<double> vel_des{vettore[NJ*campioni+time_step * NJ +0], vettore[NJ*campioni+time_step * NJ +1], vettore[NJ*campioni+time_step * NJ +2],vettore[NJ*campioni+time_step * NJ +3], vettore[NJ*campioni+time_step * NJ +4], vettore[NJ*campioni+time_step * NJ +5], vettore[NJ*campioni+time_step * NJ +6]};
				traj_msg.velocity = vel_des;
				std::vector<double> acc_des{vettore[2*NJ*campioni+time_step * NJ +0], vettore[2*NJ*campioni+time_step * NJ +1], vettore[2*NJ*campioni+time_step * NJ +2],vettore[2*NJ*campioni+time_step * NJ +3], vettore[2*NJ*campioni+time_step * NJ +4], vettore[2*NJ*campioni+time_step * NJ +5], vettore[2*NJ*campioni+time_step * NJ +6]};
				traj_msg.effort = acc_des;
				pub_cmd.publish(traj_msg);
			time_step++;

		}
		}
		
	}
	return 0;
}