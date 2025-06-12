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
#include <eigen3/Eigen/Geometry> // Per AngleAxisd, Quaterniond

// Funzione per calcolare l'errore di orientamento
Eigen::Vector3d getOrientationError(const Eigen::Quaterniond &q_desired, const Eigen::Quaterniond &q_current)
{
	Eigen::Quaterniond q_err_conj = q_current.conjugate(); // q_current.inverse() se normalizzato
	Eigen::Quaterniond q_error = q_desired * q_err_conj;
	if (q_error.w() < 0)
		q_error.coeffs() *= -1;
	q_error.normalize();

	Eigen::AngleAxisd angle_axis_error(q_error);
	  // Protezione numerica su angoli molto piccoli
    if (angle_axis_error.angle() < 1e-8)
        return Eigen::Vector3d::Zero();
	// L'errore è spesso rappresentato come angle * axis.
	// Se l'angolo è piccolo, questo approssima 2 * q_error.vec() (parte vettoriale del quaternione)
	return angle_axis_error.angle() * angle_axis_error.axis();
}

// Solutore IK iterativo (Damped Least Squares)
bool solveIK_DLS(
	thunder_franka &robot, // Passa per riferimento per poter chiamare set_q
	const Eigen::Vector3d &target_pos,
	const Eigen::Quaterniond &target_orient,
	const Eigen::Matrix<double, 7, 1> &q_initial_guess,
	Eigen::Matrix<double, 7, 1> &q_solution,
	const double q_lim_low[],			 // Passa i limiti
	const double q_lim_upp[],			 // Passa i limiti
	int max_iterations = 200,			 // Numero massimo di iterazioni
	double position_tolerance = 1e-4,	 // 0.1 mm
	double orientation_tolerance = 1e-3, // ~0.05 gradi (rad)
	double lambda_damping = 0.1,		 // Fattore di smorzamento
	double alpha_step = 0.5)			 // Dimensione del passo
{
	q_solution = q_initial_guess;
	int NJ = robot.get_numJoints(); // Ottieni il numero di giunti

	Eigen::Matrix<double, 6, 1> error_vector;
	Eigen::Matrix<double, 6, Eigen::Dynamic> J_ee(6, NJ);
	Eigen::Matrix<double, Eigen::Dynamic, 1> delta_q(NJ, 1), delta_q_giosto(NJ, 1);

	ROS_INFO("Starting IK. Initial q: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]",
			 q_solution(0), q_solution(1), q_solution(2), q_solution(3), q_solution(4), q_solution(5), q_solution(6));

	for (int iter = 0; iter < max_iterations; ++iter)
	{
		// 1. Imposta lo stato cinematico corrente nel robot per FK e Jacobiana
		robot.set_q(q_solution); // metodo setter!

		// 2. Calcola la posa corrente dell'EE (FK)
		Eigen::MatrixXd T_0_ee_mat = robot.get_T_0_ee(); // Chiama dopo set_q
		Eigen::Isometry3d current_transform = Eigen::Isometry3d::Identity();
		current_transform.matrix() = T_0_ee_mat; // Converte MatrixXd in Isometry3d

		Eigen::Vector3d current_pos = current_transform.translation();
		Eigen::Quaterniond current_orient(current_transform.rotation());
		current_orient.normalize();

		// 3. Calcola l'errore (posizione e orientamento)
		error_vector.head(3) = target_pos - current_pos;
		error_vector.tail(3) = getOrientationError(target_orient, current_orient);

		// 4. Controlla la convergenza
		if (error_vector.head(3).norm() < position_tolerance &&
			error_vector.tail(3).norm() < orientation_tolerance)
		{
			ROS_INFO("IK converged in %d iterations.", iter + 1);
			ROS_INFO("Final q: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]",
					 q_solution(0), q_solution(1), q_solution(2), q_solution(3), q_solution(4), q_solution(5), q_solution(6));
			ROS_INFO("Final Pos Error: %.5f m, Orient Error: %.5f rad", error_vector.head(3).norm(), error_vector.tail(3).norm());
			return true;
		}

		// 5. Ottieni la Jacobiana
		J_ee = robot.get_J_ee(); // Chiama dopo set_q

		// 6. Calcola il passo dei giunti (DLS: J_pinv_dls = J^T * (J * J^T + lambda^2 * I)^-1 )
		// Eigen::MatrixXd I_6x6 = Eigen::MatrixXd::Identity(6, 6);
		// Eigen::MatrixXd JJT_lambdaI = J_ee * J_ee.transpose() + lambda_damping * lambda_damping * I_6x6;
		// delta_q = J_ee.transpose() * JJT_lambdaI.inverse() * error_vector;
		delta_q = robot.get_J_ee_pinv() * error_vector;

		// 7. Aggiorna gli angoli di giunto
		q_solution += alpha_step * delta_q;

		// 8. Applica i limiti di giunto
		for (int j = 0; j < NJ; ++j)
		{
			q_solution(j) = std::max(q_lim_low[j], std::min(q_solution(j), q_lim_upp[j]));
		}
	}

	ROS_ERROR("IK failed to converge after %d iterations.", max_iterations);
	ROS_INFO("Last q: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]",
			 q_solution(0), q_solution(1), q_solution(2), q_solution(3), q_solution(4), q_solution(5), q_solution(6));
	ROS_INFO("Last Pos Error: %.5f m, Orient Error: %.5f rad", error_vector.head(3).norm(), error_vector.tail(3).norm());
	return false;
}

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
	// cout << "NJ" << NJ << endl;
	Eigen::VectorXd q(NJ), dq(NJ), dqr(NJ), ddqr(NJ);

	ros::Publisher pub_cmd = node_handle.advertise<sensor_msgs::JointState>("/computed_torque_controller/command", 1000);
	ros::Subscriber sub_joints = node_handle.subscribe<sensor_msgs::JointState>("/franka_state_controller/joint_states", 1, &jointsCallback);
	// ros::Subscriber sub_pose =  node_handle.subscribe("/franka_state_controller/franka_ee_pose", 1, &poseCallback);

	// creating trajectory message
	sensor_msgs::JointState traj_msg;

	// SET SLEEP TIME 1000 ---> 1 kHz
	double frequenza = 10; // Hz
	ros::Rate loop_rate(frequenza);
	ros::Rate loop_rate_controller(200); // Hz

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
	double tc = 0;
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
			cout << "choice:   (1: joints min-jerk,  2: go toDouble() init,  3: go toDouble() random  4: yaml 6: optimal_min-jerk 7:go to (x y z r p y) " << endl;
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

			// q_int << -1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972; // Posizioni iniziali

			for (int i = 0; i < 7; i++)
			{
				double q_low = q_lim_low[i];
				double q_upp = q_lim_upp[i];
				qf(i) = q_low + (float(rand()) / RAND_MAX) * (q_upp - q_low);
			}

			std::cout << "Random final joint positions: " << qf.transpose() << std::endl;

			// qf << -1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972; // Posizioni finali

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
		else if (choice == 7) // OPZIONE PER POSA EE
		{
			Eigen::Vector3d target_ee_pos_input;
			Eigen::Quaterniond target_ee_orient_input;
			// Per semplicità, chiediamo angoli RPY (in gradi) e li convertiamo
			double roll_deg, pitch_deg, yaw_deg;

			cout << "Enter desired EE position (x y z) in meters: ";
			cin >> target_ee_pos_input.x() >> target_ee_pos_input.y() >> target_ee_pos_input.z();

			cout << "Enter desired EE orientation RPY (roll pitch yaw) in degrees: ";
			cin >> roll_deg >> pitch_deg >> yaw_deg;

			// Converti gradi in radianti
			double roll_rad = roll_deg * M_PI / 180.0;
			double pitch_rad = pitch_deg * M_PI / 180.0;
			double yaw_rad = yaw_deg * M_PI / 180.0;

			// Converti RPY in Quaternione (convenzione ZYX per RPY)
			Eigen::AngleAxisd rollAngle(roll_rad, Eigen::Vector3d::UnitX());
			Eigen::AngleAxisd pitchAngle(pitch_rad, Eigen::Vector3d::UnitY());
			Eigen::AngleAxisd yawAngle(yaw_rad, Eigen::Vector3d::UnitZ());
			target_ee_orient_input = yawAngle * pitchAngle * rollAngle;
			target_ee_orient_input.normalize();

			cout << "Target EE Position: " << target_ee_pos_input.transpose() << endl;
			cout << "Target EE Orientation (Quaternion w,x,y,z): "
				 << target_ee_orient_input.w() << ", "
				 << target_ee_orient_input.x() << ", "
				 << target_ee_orient_input.y() << ", "
				 << target_ee_orient_input.z() << endl;

			cout << "Enter duration (tf) for the movement: ";
			cin >> tf;

			ros::spinOnce(); //  q0 (stato attuale dei giunti) aggiornato
			if (!init_q0)
			{
				ROS_ERROR("Initial joint states not received yet. Cannot solve IK.");
				continue;
			}

			Eigen::Matrix<double, 7, 1> q_target_ik;
			bool ik_solved = false;

			// Chiama il solutore IK
			// 'robot' è l'istanza della tua classe thunder_franka
			ik_solved = solveIK_DLS(robot, // L'oggetto robot
									target_ee_pos_input,
									target_ee_orient_input,
									q0,			 // Stima iniziale (configurazione corrente)
									q_target_ik, // Output: soluzione q
									q_lim_low,	 // Limiti inferiori dei giunti
									q_lim_upp);	 // Limiti superiori dei giunti

			if (ik_solved)
			{
				qf = q_target_ik; // Imposta la configurazione finale dei giunti
				ROS_INFO_STREAM("IK successful. Target joint configuration: " << qf.transpose());

				// Ora la logica di interpolazione esistente prenderà qf come target
				choice = 6; // Imposta choice a 6 per usare l'interpolazione min-jerk

				ROS_INFO_STREAM("TF for movement: " << tf << " seconds");

				v0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità iniziali

				vf << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità finali

				a0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni iniziali

				af << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni finali
			}
			else
			{
				ROS_ERROR("Failed to find IK solution for the desired pose. Skipping movement.");
				// Stampa l'ultima configurazione q0 e la posa target per debugging
				ROS_INFO_STREAM("Current q0 for IK: " << q0.transpose());
				ROS_INFO_STREAM("Target EE Pos: " << target_ee_pos_input.transpose());
				ROS_INFO_STREAM("Target EE Orient (quat w,x,y,z): " << target_ee_orient_input.w() << ", " << target_ee_orient_input.vec().transpose());
				continue; // Torna al menu
			}
		}

		if (!init_flag)
		{
			init_q0 = false;
			init_flag = true;
		}

		int campioni = tf * frequenza + 1; // Numero di campioni

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

			ros::Duration(1.0 / frequenza).sleep(); // pausa forzata all'inizio affinchè si possa calcolare la traiettoria adeguatamente

			q_int = q0;

			// cout << "q0 "<< q0 << endl;
			// cout << "q_int " << q_int << endl;

			// Calcolo dei coefficienti per ciascun giunto
			std::vector<std::vector<double>> joint_coeffs(NJ);
			for (int i = 0; i < NJ; ++i)
			{
				joint_coeffs[i] = calculateCoefficients(q_int[i], qf[i], v0[i], vf[i], a0[i], af[i], t_init.toSec(), tf);
				// std::cout << "coeffs: " << joint_coeffs[i][0] << ", " << joint_coeffs[i][1] << ", " << joint_coeffs[i][2] << ", " << joint_coeffs[i][3] << ", " << joint_coeffs[i][4] << ", " << joint_coeffs[i][5] << std::endl;
			}
			while (t <= tf)
			{
				// std::cout << "Time_inizio: " << t << std::endl;
				//  Generazione della traiettoria
				for (int i = 0; i < NJ; i++)
				{

					double pos, vel, acc;
					calculateTrajectory(t, t_init.toSec(), joint_coeffs[i], pos, vel, acc);
					// std::cout << "  Joint " << i << ": Pos = " << pos << ", Vel = " << vel << ", Acc = " << acc << std::endl;
					//   // In riga tutti i giunti e colonn i vari step
					//   POS(i, time_step) = pos;
					//   VEL(i, time_step) = vel;
					//   ACC(i, time_step) = acc;

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

				// std::vector<double> pos_des{POS_INIT[0], POS_INIT[1], POS_INIT[2], POS_INIT[3], POS_INIT[4], POS_INIT[5], POS_INIT[6]};
				// traj_msg.position = pos_des;
				// std::vector<double> vel_des{VEL_INIT[0], VEL_INIT[1], VEL_INIT[2], VEL_INIT[3], VEL_INIT[4], VEL_INIT[5], VEL_INIT[6]};
				// traj_msg.velocity = vel_des;
				// std::vector<double> acc_des{ACC_INIT[0], ACC_INIT[1], ACC_INIT[2], ACC_INIT[3], ACC_INIT[4], ACC_INIT[5], ACC_INIT[6]};
				// traj_msg.effort = acc_des;
				// pub_cmd.publish(traj_msg);

				loop_rate.sleep();

				t = (ros::Time::now() - t_init).toSec();
			}

			OptimizationData optData;
			optData.robot; // Inizializza l'oggetto robot
			optData.q0 = q_int;
			optData.qf = qf;
			optData.v0 = v0;
			optData.a0 = a0;
			optData.vf = vf;
			optData.af = af;
			optData.size_q = size_q;
			optData.dt = 1.0 / frequenza; // Passo temporale
			optData.campioni = campioni;

			// define q0,dq, ddq as 7x1 matrix
			std::vector<double> ub(3 * NJ * campioni), lb(3 * NJ * campioni), ubq(NJ), lbq(NJ), ubdq(NJ), lbdq(NJ), ubddq(NJ), lbddq(NJ);
			lbq = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
			ubq = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
			lbdq = {-2.175, -2.175, -2.175, -2.175, -2.61, -2.61, -2.61};
			ubdq = {2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61};
			lbddq = {-15, -7.5, -10, -12.5, -15, -20, -20};
			ubddq = {15, 7.5, 10, 12.5, 15, 20, 20};

			// define the optimization problem
			nlopt::opt opt(nlopt::LD_MMA, 3 * NJ * campioni);
			opt.set_min_objective(objective, &optData);
			opt.set_xtol_rel(1e-4);

			float tol = 4e-2;
			// Vincoli sulle condizioni iniziali

			for (int i = 0; i < NJ; i++)
			{
				lb[i] = q_int[i] - tol;
				ub[i] = q_int[i] + tol;
				lb[size_q + i] = v0[i] - tol;
				ub[size_q + i] = v0[i] + tol;
				lb[2 * size_q + i] = a0[i] - tol;
				ub[2 * size_q + i] = a0[i] + tol;
			}

			// Vincoli sulle condizioni intermedie

			for (int j = 1; j < campioni - 1; j++)
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

			// Vincoli sulle condizioni finali

			for (int i = 0; i < NJ; i++)
			{
				lb[(campioni - 1) * NJ + i] = qf[i] - tol;
				// std::cout << "lb " << campioni-1 * NJ + i << ":" << lb[campioni-1 * NJ + i] << std::endl;
				ub[(campioni - 1) * NJ + i] = qf[i] + tol;
				// std::cout << "ub " << campioni-1 * NJ + i << ":" << ub[campioni-1 * NJ + i] << std::endl;
				lb[(campioni - 1) * NJ + size_q + i] = vf[i] - tol;
				// std::cout << "lb " << campioni-1 * NJ + size_q + i << ":" << lb[campioni-1 * NJ + size_q + i] << std::endl;
				ub[(campioni - 1) * NJ + size_q + i] = vf[i] + tol;
				// std::cout << "ub " << campioni-1 * NJ + size_q + i << ":" << ub[campioni-1 * NJ + size_q + i] << std::endl;
				lb[(campioni - 1) * NJ + 2 * size_q + i] = af[i] - tol;
				// std::cout << "lb " << campioni-1 * NJ + 2 * size_q + i << ":" << lb[campioni-1 * NJ + 2 * size_q + i] << std::endl;
				ub[(campioni - 1) * NJ + 2 * size_q + i] = af[i] + tol;
				// std::cout << "ub " << campioni-1 * NJ + 2 * size_q + i << ":" << ub[campioni-1 * NJ + 2 * size_q + i] << std::endl;
			}

			opt.set_upper_bounds(ub);
			opt.set_lower_bounds(lb);

			std::vector<ConsistencyConstraintIneq> constraints;
			std::vector<std::shared_ptr<ObstacleConstraintIneq>> sphere_constraints;
			const double eps = 4e-2;

			int numero_totale_vincoli = (campioni - 1) * 7; // <-- Calcola il numero totale

			const double r_s = 0.05;	   // raggio ostacolo
			const double d_safe = 0.1; // margine sicurezza
			Eigen::Vector3d p_ostacolo(0.25, 0.25, 0.8);

			if (numero_totale_vincoli > 0)
			{ // evita di chiamare reserve(0)
				constraints.reserve(numero_totale_vincoli);
			}

			for (int k = 0; k < campioni - 1; k++)
			{
				// Posizione
				constraints.push_back({k, NJ, size_q, optData.dt, 0, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				constraints.push_back({k, NJ, size_q, optData.dt, 0, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				// Velocità
				constraints.push_back({k, NJ, size_q, optData.dt, 1, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				constraints.push_back({k, NJ, size_q, optData.dt, 1, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				// Accelerazione
				constraints.push_back({k, NJ, size_q, optData.dt, 2, +1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				constraints.push_back({k, NJ, size_q, optData.dt, 2, -1});
				opt.add_inequality_constraint(consistency_ineq, &constraints.back(), eps);

				// --- Vincolo di evitamento ostacolo (sfera) ---
				auto c = std::make_shared<ObstacleConstraintIneq>(k, NJ, r_s, d_safe, p_ostacolo, robot);
				sphere_constraints.push_back(c);
				opt.add_inequality_constraint(avoid_sphere, c.get(), eps);

				// sphere_constraints.push_back({k, NJ, r_s, d_safe, p_ostacolo, robot});
				// opt.add_inequality_constraint(avoid_sphere, &sphere_constraints.back(), eps);
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

			// Traiettoria fianale smussata
			int count = 1;
			for (int j = 0; j < campioni - 1; j++)
			{
				double t_start = ros::Time::now().toSec();
				// Calcolo dei coefficienti per ciascun giunto
				std::vector<std::vector<double>> joint_coeffs(NJ);
				for (int i = 0; i < NJ; ++i)
				{
					q_int(i) = vettore[j * NJ + i];
					qf(i) = vettore[(j + 1) * NJ + i];
					v0(i) = vettore[j * NJ + size_q + i];
					vf(i) = vettore[(j + 1) * NJ + size_q + i];
					a0(i) = vettore[j * NJ + 2 * size_q + i];
					af(i) = vettore[(j + 1) * NJ + 2 * size_q + i];
					tf = t_start + 1.0 / frequenza;
					joint_coeffs[i] = calculateCoefficients(q_int[i], qf[i], v0[i], vf[i], a0[i], af[i], t_start, tf);
					// std::cout << "coeffs: " << joint_coeffs[i][0] << ", " << joint_coeffs[i][1] << ", " << joint_coeffs[i][2] << ", " << joint_coeffs[i][3] << ", " << joint_coeffs[i][4] << ", " << joint_coeffs[i][5] << std::endl;
				}

				t = t_start;
				// std::cout << "tempo_inizio" << t_start << std::endl;
				// std::cout << "tempo_fine" << tf << std::endl;

				while (t <= tf)
				{

					for (int i = 0; i < NJ; i++)
					{

						double pos, vel, acc;
						calculateTrajectory(t, t_init.toSec(), joint_coeffs[i], pos, vel, acc);
						// std::cout << "  Joint " << i << ": Pos = " << pos << ", Vel = " << vel << ", Acc = " << acc << std::endl;
						//  // In riga tutti i giunti e colonn i vari step
						//  POS(i, time_step) = pos;
						//  VEL(i, time_step) = vel;
						//  ACC(i, time_step) = acc;

						// Incolonna e colleziona tutte le varibili di giunto per ogni step
						POS_INIT(j * NJ + i) = pos;
						VEL_INIT(j * NJ + i) = vel;
						ACC_INIT(j * NJ + i) = acc;

						// // Incolonna tutte le varibili di giunto ad ogni step
						// POS_INIT(i) = pos;
						// VEL_INIT(i) = vel;
						// ACC_INIT(i) = acc;
					}
					std::vector<double> pos_des{POS_INIT[j * NJ + 0], POS_INIT[j * NJ + 1], POS_INIT[j * NJ + 2], POS_INIT[j * NJ + 3], POS_INIT[j * NJ + 4], POS_INIT[j * NJ + 5], POS_INIT[j * NJ + 6]};
					traj_msg.position = pos_des;
					std::vector<double> vel_des{VEL_INIT[j * NJ + 0], VEL_INIT[j * NJ + 1], VEL_INIT[j * NJ + 2], VEL_INIT[j * NJ + 3], VEL_INIT[j * NJ + 4], VEL_INIT[j * NJ + 5], VEL_INIT[j * NJ + 6]};
					traj_msg.velocity = vel_des;
					std::vector<double> acc_des{ACC_INIT[j * NJ + 0], ACC_INIT[j * NJ + 1], ACC_INIT[j * NJ + 2], ACC_INIT[j * NJ + 3], ACC_INIT[j * NJ + 4], ACC_INIT[j * NJ + 5], ACC_INIT[j * NJ + 6]};
					traj_msg.effort = acc_des;
					pub_cmd.publish(traj_msg);

					loop_rate_controller.sleep();

					t = ros::Time::now().toSec();
					// std::cout << "Time: " << t << std::endl;

					if (t > tf)
					{
						// std::cout << "fine cicl:  " << count << std::endl;
						count++;
					}
				}
			}

			// // Print the optimized values
			// printf("Optimized values:\n");
			// printf("q_opt:\n");
			// for (int i = 0; i < POS_INIT.size(); i++)
			// {
			// 	printf("%g ", vettore[i]);
			// 	printf(",\n");
			// }
			// printf(",\n");
			// printf("dq_opt:\n");
			// for (int i = POS_INIT.size(); i < 2 * POS_INIT.size(); i++)
			// {
			// 	printf("%g ", vettore[i]);
			// 	printf(",\n");
			// }
			// printf(",\n");
			// printf("ddq_opt:\n");
			// for (int i = 2 * POS_INIT.size(); i < 3 * POS_INIT.size(); i++)
			// {
			// 	printf("%g ", vettore[i]);
			// 	printf(",\n");
			// }

			// // time_step = 0;

			// // while (time_step < campioni)

			// // {

			// // 	std::vector<double> pos_des{vettore[time_step * NJ + 0], vettore[time_step * NJ + 1], vettore[time_step * NJ + 2], vettore[time_step * NJ + 3], vettore[time_step * NJ + 4], vettore[time_step * NJ + 5], vettore[time_step * NJ + 6]};

			// // 	traj_msg.position = pos_des;
			// // 	std::vector<double> vel_des{vettore[NJ * campioni + time_step * NJ + 0], vettore[NJ * campioni + time_step * NJ + 1], vettore[NJ * campioni + time_step * NJ + 2], vettore[NJ * campioni + time_step * NJ + 3], vettore[NJ * campioni + time_step * NJ + 4], vettore[NJ * campioni + time_step * NJ + 5], vettore[NJ * campioni + time_step * NJ + 6]};
			// // 	traj_msg.velocity = vel_des;
			// // 	std::vector<double> acc_des{vettore[2 * NJ * campioni + time_step * NJ + 0], vettore[2 * NJ * campioni + time_step * NJ + 1], vettore[2 * NJ * campioni + time_step * NJ + 2], vettore[2 * NJ * campioni + time_step * NJ + 3], vettore[2 * NJ * campioni + time_step * NJ + 4], vettore[2 * NJ * campioni + time_step * NJ + 5], vettore[2 * NJ * campioni + time_step * NJ + 6]};
			// // 	traj_msg.effort = acc_des;
			// // 	pub_cmd.publish(traj_msg);

			// // 	time_step++;
			// // }
		}
	}

	return 0;
}