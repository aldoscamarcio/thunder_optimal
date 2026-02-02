#include <iostream>
#include <eigen3/Eigen/Dense>
#include <unistd.h>
#include <cstdlib>
#include <signal.h>
#include <math.h>
#include <vector>
#include <geometry_msgs/PoseStamped.h>
#include "ros/ros.h"
#include <random>
#include "utils/thunder_franka.h"
#include "utils/thunder_optimization.h"
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <chrono>
#include <iostream>
#include <iomanip>
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
#include "collision/CollisionEngine.hpp"
#include <rosbag/bag.h>
#include <std_msgs/Float64.h>

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
	int max_iterations = 1000,			 // Numero massimo di iterazioni
	double position_tolerance = 1e-3,	 // 0.1 mm
	double orientation_tolerance = 1e-2, // ~0.05 gradi (rad)
	double lambda_damping = 0.5,		 // Fattore di smorzamento
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

std::string nloptResultToString(nlopt::result res) {
    switch (res) {
        case nlopt::SUCCESS: return "SUCCESS (Convergenza raggiunta)";
        case nlopt::STOPVAL_REACHED: return "STOPVAL_REACHED";
        case nlopt::FTOL_REACHED: return "FTOL_REACHED (Variazione f minima)";
        case nlopt::XTOL_REACHED: return "XTOL_REACHED (Variazione x minima)";
        case nlopt::MAXEVAL_REACHED: return "MAXEVAL_REACHED (Troppe iterazioni)";
        case nlopt::MAXTIME_REACHED: return "MAXTIME_REACHED (Tempo scaduto)";
        case nlopt::FAILURE: return "FAILURE (Errore generico)";
        case nlopt::INVALID_ARGS: return "INVALID_ARGS (Argomenti errati)";
        case nlopt::OUT_OF_MEMORY: return "OUT_OF_MEMORY";
        case nlopt::ROUNDOFF_LIMITED: return "ROUNDOFF_LIMITED (Errori arrotondamento)";
        case nlopt::FORCED_STOP: return "FORCED_STOP";
        default: return "UNKNOWN CODE";
    }
}

Eigen::Vector3d generaOstacoloRandom() {
	// Inizializza generatore casuale con seed fisso per riproducibilità
	std::mt19937 gen(1);
	
	// Definisci range per le coordinate (modifica secondo necessità)
	std::uniform_real_distribution<double> dist_x(-0.7, 0.7);   // x tra 0.0 e 0.2
	std::uniform_real_distribution<double> dist_y(-0.7, 0.7); // y tra -0.5 e -0.2
	std::uniform_real_distribution<double> dist_z(0.06, 1.2);   // z tra 0.4 e 0.6
	
	return Eigen::Vector3d(dist_x(gen), dist_y(gen), dist_z(gen));
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

int main(int argc, char **argv)
{
	ros::init(argc, argv, "menu");

	ros::NodeHandle node_handle;
	thunder_franka robot;
	robot.load_conf(conf_file);
	int NJ = robot.get_numJoints();
	Eigen::VectorXd q(NJ), dq(NJ), dqr(NJ), ddqr(NJ);

	ros::Publisher pub_cmd = node_handle.advertise<sensor_msgs::JointState>("/computed_torque_controller/command", 1000);
	ros::Subscriber sub_joints = node_handle.subscribe<sensor_msgs::JointState>("/franka_state_controller/joint_states", 1, &jointsCallback);
	ros::Publisher capsule_viz_pub = node_handle.advertise<visualization_msgs::MarkerArray>("robot_capsules_viz", 10);
	ros::Publisher marker_pub = node_handle.advertise<visualization_msgs::MarkerArray>("/optimization_markers", 10);
	// ros::Publisher path_pub = node_handle.advertise<nav_msgs::Path>("/end_effector_path", 1);
	// ros::Subscriber sub_pose =  node_handle.subscribe("/franka_state_controller/franka_ee_pose", 1, &poseCallback);

	// creating trajectory message
	sensor_msgs::JointState traj_msg;

	// SET SLEEP TIME 1000 ---> 1 kHz
	double frequenza = 10;
	ros::Rate loop_rate(frequenza);
	double frequenza_controller = 200;
	ros::Rate loop_rate_controller(frequenza_controller);

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

	// ===========================================
	// CAPSULE GENERATE DA fr3_franka_hand.urdf
	// ===========================================
	std::vector<Capsule> capsule_definitions;
	{
		Capsule cap;
		cap.link_index = 0;
		cap.radius = 0.055000;
		cap.length = 0.030000;
		cap.T_offset << 0.0000, 0.0000, 1.0000, -0.0750,
			0.0000, 1.0000, 0.0000, 0.0000,
			-1.0000, 0.0000, 0.0000, 0.0600,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link1 (Index 1) ---
	{
		Capsule cap;
		cap.link_index = 1;
		cap.radius = 0.060000;
		cap.length = 0.183000;
		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 1.0000, -0.1915,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link2 (Index 2) ---
	{
		Capsule cap;
		cap.link_index = 2;
		cap.radius = 0.070000;
		cap.length = 0.240000;
		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 1.0000, 0.0000,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link3 (Index 3) ---
	{
		Capsule cap;
		cap.link_index = 3;
		cap.radius = 0.090000;
		cap.length = 0.150000;
		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 1.0000, -0.1450,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link4 (Index 4) ---
	{
		Capsule cap;
		cap.link_index = 4;
		cap.radius = 0.080000;
		cap.length = 0.200000;
		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 1.0000, 0.0000,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link5 (Index 5) ---
	{
		Capsule cap;
		cap.link_index = 5;
		cap.radius = 0.090000;
		cap.length = 0.140000;
		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 1.0000, -0.2600,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	{
		Capsule cap;
		cap.link_index = 5;
		cap.radius = 0.055000;
		cap.length = 0.150000;
		cap.T_offset << 0.9968, -0.0799, 0.0000, 0.0000,
			0.0799, 0.9968, 0.0000, 0.0800,
			0.0000, 0.0000, 1.0000, -0.1300,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link6 (Index 6) ---
	{
		Capsule cap;
		cap.link_index = 6;
		cap.radius = 0.070000;
		cap.length = 0.150000;
		// Offset X corretto a -0.0100 per centrare sul giunto
cap.T_offset << 1.0000,  0.0000,  0.0000, -0.0100,
                0.0000,  0.0000, -1.0000,  0.0300,  // Y → -Z
                0.0000,  1.0000,  0.0000, -0.0100,  // Z → Y
                0.0000,  0.0000,  0.0000,  1.0000;
		capsule_definitions.push_back(cap);
	}

	// --- fr3_link7 (Index 7)
	{
		// Parte 1: Braccio Orizzontale
		Capsule cap;
		cap.link_index = 7;
		cap.radius = 0.050000;
		cap.length = 0.120000;

		cap.T_offset << 0.0000, 0.0000, 1.0000, 0.0600,
			0.0000, 1.0000, 0.0000, 0.0000,
			-1.0000, 0.0000, 0.0000, 0.0100,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

	{
		// Parte 2: Flangia Verticale/Obliqua
		Capsule cap;
		cap.link_index = 7;
		cap.radius = 0.060000;
		cap.length = 0.190000;

		cap.T_offset << 1.0000, 0.0000, 0.0000, 0.1050,
			0.0000, 0.0000, -1.0000, -0.0100,
			0.0000, 1.0000, 0.0000, 0.0000,
			0.0000, 0.0000, 0.0000, 1.0000;
		capsule_definitions.push_back(cap);
	}

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
			cout << "choice:   (1: joints min-jerk,  2: go to init,  3: go to random  4: yaml 6: optimal_min-jerk 7: go to (x y z r p y)" << endl;
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
			tf = 5.0;
			for (int i = 0; i < 7; i++)
			{
				double q_low = q_lim_low[i];
				double q_upp = q_lim_upp[i];
				qf(i) = q_low + (float(rand()) / RAND_MAX) * (q_upp - q_low);
			}

			std::cout << "Random final joint positions: " << qf.transpose() << std::endl;
			v0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità iniziali
			vf << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Velocità finali
			a0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni iniziali
			af << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // Accelerazioni finali
		}
		else if (choice == 7) // Opzione per Posizione e Orientamento EE
		{
			Eigen::Vector3d target_ee_pos_input;
			Eigen::Quaterniond target_ee_orient_input;
			double roll_deg, pitch_deg, yaw_deg;

			// target_ee_pos_input.x() = 0.30; // Posizione EE desiderata in metri
			// target_ee_pos_input.y() = 0.40;
			// target_ee_pos_input.z() = 0.68; // Posizione EE desiderata in metri

			// RANDOM Z per test Autocollisione
			// target_ee_pos_input.z() = -0.3 + (float(rand()) / RAND_MAX) * (0.8 - (-0.3)); // Tra -0.3 e 0.8 m

			roll_deg = -180; // Angolo roll in gradi
			pitch_deg = 0;   // Angolo pitch in gradi
			yaw_deg = 60;	   // Angolo yaw in gradi

			// Sezione di input manuale
			// cout << "Enter desired EE position (x y z) in meters: ";
			// cin >> target_ee_pos_input.x() >> target_ee_pos_input.y() >> target_ee_pos_input.z();

			// cout << "Enter desired EE orientation RPY (roll pitch yaw) in degrees: ";
			// cin >> roll_deg >> pitch_deg >> yaw_deg;

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

			cout << "Do you want to optimize the movement? (1: Yes, 0: No): ";
			int optimize_movement;
			cin >> optimize_movement;

			// Chiedi durata del movimento
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
			ik_solved = solveIK_DLS(robot, // L'oggetto robot
									target_ee_pos_input,
									target_ee_orient_input,
									q0,			 // Stima iniziale (configurazione corrente)
									q_target_ik, // Output: soluzione q
									q_lim_low,	 // Limiti inferiori dei giunti
									q_lim_upp);	 // Limiti superiori dei giunti

			// 1. Assegnazione Target (Comune a entrambi i casi)
			// qf = ik_solved ? q_target_ik : qf; // Scommenta se usi la IK calcolata quando disponibile
			// qf << 0.48145, 0.482112, 0.776829, -1.50611, -0.254019, 1.85678, 0.178126; // Target MultiOstacolo
			// qf << 0.5889, 0.3757, 0.9471, -1.6643, -0.2357, 1.9641, 0.1780;
			qf << 0.335, -0.321, 0.701, -2.145, 0.174, 1.884, 1.740; // Posizione di Arrivo test singolo ostacolo

			if (ik_solved)
			{
				ROS_INFO_STREAM("IK successful. Target joint configuration: " << qf.transpose());
			}
			else
			{
				ROS_ERROR("Failed to find IK solution. Forcing movement anyway.");
				// qf = q_target_ik;
				qf << 0.335, -0.321, 0.701, -2.145, 0.174, 1.884, 1.740; // Posizione di Arrivo test Collisione
				// qf << 0.48145, 0.482112, 0.776829, -1.50611, -0.254019, 1.85678, 0.178126; // Posizione di Arrivo test MultiOstacolo 
				// qf << 0.5889, 0.3757, 0.9471, -1.6643, -0.2357, 1.9641, 0.1780;
				ROS_INFO_STREAM("Debug - Current q0: " << q0.transpose());
				ROS_INFO_STREAM("Debug - Target EE: " << target_ee_pos_input.transpose());
				ROS_INFO_STREAM("Force Target: " << qf.transpose());
			}

			if (optimize_movement == 1)
			{
				choice = 6; // Min-Jerk Optimization
				// Reset dinamica (usa setZero di Eigen per pulizia)
				v0.setZero();
				vf.setZero();
				a0.setZero();
				af.setZero();
			}
			else if (optimize_movement == 0)
			{
				choice = 1; // Interpolazione semplice
			}
			else
			{
				continue;
			}
		}

		if (!init_flag)
		{
			init_q0 = false;
			init_flag = true;
		}

		int campioni = tf * frequenza + 1; // Numero di campioni

		int size_q = NJ * campioni; // Dimensione di q

		Eigen::VectorXd POS_INIT(NJ * campioni), VEL_INIT(NJ * campioni), ACC_INIT(NJ * campioni), VEL(NJ * campioni), POS(NJ * campioni), ACC(NJ * campioni);

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

			// 1. Prepara i dati
			JointLimitsData pos_data;
			pos_data.NJ = NJ;
			pos_data.dt = optData.dt;    // Tuo dt
			pos_data.q0 = q0; // Start conf
			pos_data.v0 = v0; // Start vel
			pos_data.qf = qf; // End acc
			pos_data.vf = vf; // End vel
			pos_data.campioni = campioni;
			pos_data.ubq = (Eigen::Matrix<double, 7, 1>() << 2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973).finished(); // Vettori limiti
			pos_data.lbq = (Eigen::Matrix<double, 7, 1>() << -2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973).finished(); // Vettori limiti
			pos_data.ubdq = (Eigen::Matrix<double, 7, 1>() << 2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61).finished();
			pos_data.lbdq = (Eigen::Matrix<double, 7, 1>() << -2.175, -2.175, -2.175, -2.175, -2.61, -2.61, -2.61).finished();

			// Imposto vettori limiti superiori e inferiori per i vincoli
			std::vector<double> ub(NJ * campioni), lb(NJ * campioni), ubq(NJ), lbq(NJ), ubdq(NJ), lbdq(NJ), ubddq(NJ), lbddq(NJ);
			lbq = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
			ubq = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
			lbdq = {-2.175, -2.175, -2.175, -2.175, -2.61, -2.61, -2.61};
			ubdq = {2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61};
			lbddq = {-15, -7.5, -10, -12.5, -15, -20, -20};
			ubddq = {15, 7.5, 10, 12.5, 15, 20, 20};

			// define the optimization problem
			nlopt::opt opt(nlopt::LD_MMA, NJ * campioni);
			opt.set_min_objective(objective, &optData);

			opt.set_ftol_rel(1e-4);	   // Tolleranza di convergenza funzione obiettivo
			// opt.set_xtol_abs(1e-5);	   // Tolleranza di convergenza variabili	Funziona bene solo 0.01 di errore
			// opt.set_xtol_rel(1e-3);	   // Tolleranza di convergenza variabili

			opt.set_param("verbosity", 1); // Verbose output
			const double tol = 1e-2;		   // Tolleranza per i vincoli di uguaglianza

			// Vincoli sulle condizioni iniziali, imponiamo all'ottimizzatore che le condizioni iniziali siano rispettate
			for (int i = 0; i < NJ; i++)
			{

				lb[i] = a0[i] - tol;
				ub[i] = a0[i] + tol;
			}

			// Vincoli sulle condizioni intermedie, solo limiti di giunto

			for (int j = 1; j < campioni - 1; j++)
			{
				for (int i = 0; i < NJ; i++)
				{

					lb[j * NJ + i] = lbddq[i];
					ub[j * NJ + i] = ubddq[i];
				}
			}

			// Vincoli sulle condizioni finali, imponiamo all'ottimizzatore che le condizioni finali siano rispettate

			for (int i = 0; i < NJ; i++)
			{
				lb[(campioni - 1) * NJ + i] = af[i] - tol;
				// std::cout << "lb " << campioni-1 * NJ + 2 * size_q + i << ":" << lb[campioni-1 * NJ + 2 * size_q + i] << std::endl;
				ub[(campioni - 1) * NJ + i] = af[i] + tol;
				// std::cout << "ub " << campioni-1 * NJ + 2 * size_q + i << ":" << ub[campioni-1 * NJ + 2 * size_q + i] << std::endl;
			}

			opt.set_upper_bounds(ub);
			opt.set_lower_bounds(lb);

			// Vincoli di consistenza + evitamento ostacolo
			std::vector<std::shared_ptr<ObstacleConstraintIneq>> sphere_constraints, plane_constraints, rectangle_constraints;
			std::vector<std::shared_ptr<SelfCollisionConstraint>> self_coll_constraints;
			// Tolleranze per i vincoli
			const double eps_sphere = 1e-4; // Tolleranza per i vincoli di evitamento ostacolo

			// Calcola il numero totale di vincoli correttamente
			int vincoli_ostacolo_per_step = 4; // Sfera + piano + self-collision + rettangolo
			int vincoli_pos_per_step = 2 * NJ; // Upper + lower per ogni giunto
			int vincoli_vel_per_step = 2 * NJ; // Upper + lower per ogni giunto
			int vincoli_per_step = vincoli_ostacolo_per_step + vincoli_pos_per_step + vincoli_vel_per_step;

			unsigned m_pos_vel = 2 * NJ * campioni;
			std::vector<double> tolm(m_pos_vel, 1e-2); // Tolleranza

			unsigned m_final =  NJ; 
			double precision_pos = 1e-4; 
			std::vector<double> tol_final(m_final, precision_pos);


			int numero_totale_vincoli = (campioni) * vincoli_per_step;
			const double r_s = 0.07;   // raggio ostacolo
			const double d_safe = 0.05; // Margine di Sicurezza (8 cm per la sfera)/*(5cm per il piano)

			Eigen::Vector3d p_sfera(0.40, 0, 0.65); // Posizione fissa ostacolo
			Eigen::Vector3d p_piano(0.0, 0.0, 0.0); // Punto sul piano
			Eigen::Vector3d p_rettangolo(0.25, -0.01, 0.50);

			Obstacle obs_sphere;
			obs_sphere.type = ObstacleType::CAPSULE;
			obs_sphere.capsule.A = p_sfera; // Centro sfera
			obs_sphere.capsule.B = p_sfera; // Stesso punto
			obs_sphere.capsule.radius = r_s;   // Raggio sfera

			Eigen::Vector3d p_sphere_blue(0.30, 0.45, 0.80);
			Obstacle obs_sphere_blue;
			const double r_sb = 0.08;   // raggio ostacolo
			const double d_safe_blue = 0.06; // Margine di Sicurezza (5 cm per la sfera blu)
    		obs_sphere_blue.type = ObstacleType::CAPSULE;
    
			obs_sphere_blue.capsule.A = p_sphere_blue;
			obs_sphere_blue.capsule.B = p_sphere_blue;
			obs_sphere_blue.capsule.radius = r_sb;   // Raggio sfera

				
			// Cilindro Verticale
			// Obstacle obs_cylinder;
			// obs_cylinder.type = ObstacleType::CAPSULE;
			// double altezza_cilindrica = 0.18;  // 10 cm tra i centri delle semisfere
			// Eigen::Vector3d p_cilindro_A(0.38, 0, 0.60 - altezza_cilindrica/2);  // (0.38, 0, 0.55)
			// Eigen::Vector3d p_cilindro_B(0.38, 0, 0.60 + altezza_cilindrica/2);  // (0.38, 0, 0.65)

			// obs_cylinder.capsule.A = p_cilindro_A;
			// obs_cylinder.capsule.B = p_cilindro_B;
			// obs_cylinder.capsule.radius = 0.05;

			// Cilindro Orizzontale
			// Per il cilindro orizzontale lungo X (capsula con estremi separati)
			// double r_c = 0.02;
			// double lunghezza_cilindro = 0.10;  // Lunghezza del cilindro (senza le semisfere)
			// double distanza_AB = lunghezza_cilindro - 2*r_c;  // Parte cilindrica pura

			// Eigen::Vector3d p_cilindro_A(0.38 - distanza_AB/2, 0, 0.60);  // Estremo sinistro
			// Eigen::Vector3d p_cilindro_B(0.38 + distanza_AB/2, 0, 0.60);  // Estremo destro

			// Obstacle obs_cilindro;
			// obs_cilindro.type = ObstacleType::CAPSULE;
			// obs_cilindro.capsule.A = p_cilindro_A;
			// obs_cilindro.capsule.B = p_cilindro_B;
			// obs_cilindro.capsule.radius = r_c;

			Obstacle obs_plane;
			obs_plane.type = ObstacleType::PLANE;
			obs_plane.plane.P0 = p_piano;				  // Punto sul piano
			obs_plane.plane.n = Eigen::Vector3d(0, 0, 1); // Normale del piano

			Obstacle obs_rectangle;
			obs_rectangle.type = ObstacleType::RECTANGLE;
			obs_rectangle.rect.P0 = p_rettangolo;												// Vertice in basso a sinistra
			obs_rectangle.rect.Ux = Eigen::Vector3d(1, 0, 0);								// Vettore direzione u (larghezza)
			obs_rectangle.rect.Uy = Eigen::Vector3d(0, 0, 1);								// Vettore direzione v (altezza)
			obs_rectangle.rect.width = 0.2;													// Lunghezza lungo u
			obs_rectangle.rect.height = 0.2;												// Lunghezza lungo v
			obs_rectangle.rect.normal = obs_rectangle.rect.Ux.cross(obs_rectangle.rect.Uy); // Normale del rettangolo

			std::vector<std::pair<int, int>> collision_pairs = {
			    {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8},
			    {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {0, 9},
			    {2, 5}, {2, 6}, {2, 7}, {2, 8}, {2, 9},
			    {3, 6}, {3, 7}, {3, 8},
			    {4, 7}, {4, 8},
				// {5, 8}, {4, 6}, {2, 4}, {3, 5}
			};

			// Preallocazione vettori vincoli

				// constraints_pos_f.clear();
				// constraints_pos_f.reserve((campioni) * vincoli_per_step + 2 * NJ);
				// constraints_vel_f.clear();
				// constraints_vel_f.reserve((campioni) * vincoli_per_step + 2 * NJ);
				// constraints_pos.clear();
				// constraints_pos.reserve(NJ * 4);
				sphere_constraints.clear();
				sphere_constraints.reserve(numero_totale_vincoli);
				plane_constraints.clear();
				plane_constraints.reserve(numero_totale_vincoli);
				self_coll_constraints.clear();
				self_coll_constraints.reserve(collision_pairs.size() * numero_totale_vincoli);
				rectangle_constraints.clear();
				rectangle_constraints.reserve(numero_totale_vincoli);

			opt.add_inequality_mconstraint(position_limits_mconstraint, &pos_data, tolm);
			opt.add_inequality_mconstraint(velocity_limits_mconstraint, &pos_data, tolm);
			opt.add_inequality_mconstraint(final_position_inequality_mconstraint, &pos_data, tol_final);
			opt.add_inequality_mconstraint(final_velocity_inequality_mconstraint, &pos_data, tol_final);


			//  Vincoli aggiuntivi per collisione con ostacoli e self-collision
			for (int k = 0; k < campioni; k++)
			{
				auto c_sphere = std::make_shared<ObstacleConstraintIneq>(
					k, NJ, d_safe, obs_sphere, &robot, optData.q0, optData.v0, optData.dt, &capsule_viz_pub);
				c_sphere->capsules_definitions = capsule_definitions;
				sphere_constraints.push_back(c_sphere);
				
				// Sfera Blu Aggiuntiva per Scenari complessi
				auto c_sphere_blue = std::make_shared<ObstacleConstraintIneq>(
					k, NJ, d_safe_blue, obs_sphere_blue, &robot, optData.q0, optData.v0, optData.dt, &capsule_viz_pub);
				c_sphere_blue->capsules_definitions = capsule_definitions;
				sphere_constraints.push_back(c_sphere_blue);

				auto c_plane = std::make_shared<ObstacleConstraintIneq>(
					k, NJ, d_safe, obs_plane, &robot, optData.q0, optData.v0, optData.dt, &capsule_viz_pub);
				c_plane->capsules_definitions = capsule_definitions;
				plane_constraints.push_back(c_plane);

				auto c_self = std::make_shared<SelfCollisionConstraint>(SelfCollisionConstraint{
					.k = k,
					.NJ = NJ,
					.dt = optData.dt,
					.d_safe = 0.10, // Margine di sicurezza per self-collision
					.q0 = optData.q0,
					.dq0 = optData.v0,
					.robot = &robot,
					.capsules_definitions = capsule_definitions,
					.collision_pairs = collision_pairs,
					.marker_pub = &marker_pub
				});
				self_coll_constraints.push_back(c_self);

				auto c_rectangle = std::make_shared<ObstacleConstraintIneq>(
					k, NJ, d_safe, obs_rectangle, &robot, optData.q0, optData.v0, optData.dt, &capsule_viz_pub);
				c_rectangle->capsules_definitions = capsule_definitions;
				rectangle_constraints.push_back(c_rectangle);
				opt.add_inequality_constraint(avoid_obstacle_generic, c_rectangle.get(), eps_sphere);

				
				// opt.add_inequality_constraint(avoid_obstacle_generic, c_sphere.get(), eps_sphere);
				// opt.add_inequality_constraint(avoid_obstacle_generic, c_sphere_blue.get(), eps_sphere);
				// opt.add_inequality_constraint(avoid_obstacle_generic, c_plane.get(), 1e-4);
				// opt.add_inequality_constraint(avoid_self_collision, c_self.get(), 1e-2);

			}

			// Define start time for optimization
			auto start_time = std::chrono::high_resolution_clock::now();

			// Define the initial guess
			std::vector<double> vettore(NJ * campioni); // Inizializza il vettore x
			for (int i = 0; i < ACC_INIT.size(); i++)
			{
				vettore[i] = ACC_INIT[i];
			}
			double minf;

			try
            {
                nlopt::result result = opt.optimize(vettore, minf);
                std::cout << "Ottimizzazione Completata con codice: " << result << std::endl;
				std::cout << nloptResultToString(result) << std::endl;
                std::cout << "Costo minimo: " << minf << std::endl;

            }
			catch (nlopt::forced_stop &e) {
				ROS_WARN("Ottimizzazione interrotta manualmente dall'utente! ");
				ROS_WARN("Salvataggio della migliore traiettoria trovata finora (Costo: %f)", minf);
				
			}
            catch (std::exception &e)
            {
                std::cerr << "Errore NLOPT!: " << e.what() << std::endl;
				continue; // Torna al menu principale in caso di errore
            }

			// Define end time for optimization
			double optimization_time = 0.0;
    		auto end_time = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::ratio<1>> elapsed_ms = end_time - start_time;

			std::cout << "  Tempo Ottimizzazione: " << std::fixed << std::setprecision(2) << elapsed_ms.count() << " s" << std::endl;
			ros::NodeHandle nh;
			// nh.setParam("/thunder/last_optimization_time", elapsed_ms.count());

			ros::param::set("/thunder/last_planner_id", "NLOPT_MMA");
        	ros::param::set("/thunder/last_planning_time_cpu", elapsed_ms.count());
        	ros::param::set("/thunder/last_optimization_time", elapsed_ms.count());

            // 1. RICOSTRUZIONE DELLA TRAIETTORIA (SENZA ESEGUIRE)
            // Inizializzione posizioni, velocità e accelerazioni
            for (int j = 0; j < NJ; ++j)
            {
                POS[j] = q0[j]; // posizione iniziale
                VEL[j] = v0[j]; // velocità iniziale
            }

            for (int k = 0; k < campioni - 1; k++)
            {
                for (int j = 0; j < NJ; ++j)
                {
                    // Qui usiamo 'vettore', quindi non doveva essere cancellato prima!
                    ACC[k * NJ + j] = vettore[k * NJ + j]; 
                    VEL[(k + 1) * NJ + j] = VEL[k * NJ + j] + ACC[k * NJ + j] * optData.dt;
                    POS[(k + 1) * NJ + j] = POS[k * NJ + j] + VEL[k * NJ + j] * optData.dt + 0.5 * ACC[k * NJ + j] * std::pow(optData.dt, 2);
                }
            }
            
            // Ora che abbiamo estratto i dati, puliamo il vettore dell'ottimizzatore
            vettore.clear();
            std::cout << "Eseguire il movimento sul robot? (y/n): ";
            
            char user_input;
            std::cin >> user_input;
            
            // Pulisce il buffer nel caso l'utente abbia premuto invio o digitato caratteri extra
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            if (user_input != 'y' && user_input != 'Y') 
            {
                std::cout << ">>> Esecuzione ANNULLATA dall'utente. Robot fermo. " << std::endl;
				continue;
            }
            else 
            {
                std::cout << "Esecuzione Traiettoria in corso" << std::endl;
                
                for (int j = 0; j < campioni - 1; j++)
                {
                    double t_start = ros::Time::now().toSec();
                    // Calcolo dei coefficienti per ciascun giunto
                    std::vector<std::vector<double>> joint_coeffs(NJ);
                    for (int i = 0; i < NJ; ++i)
                    {
                    q_int(i) = POS[j * NJ + i];
                    qf(i) = POS[(j + 1) * NJ + i];
                    v0(i) = VEL[j * NJ + i];
                    vf(i) = VEL[(j + 1) * NJ + i];
                    a0(i) = ACC[j * NJ + i];
                    af(i) = ACC[(j + 1) * NJ + i];
                    tf = t_start + 1.0 / frequenza;

                    joint_coeffs[i] = calculateCoefficients(q_int[i], qf[i], v0[i], vf[i], a0[i], af[i], t_start, tf);
                    }

                    t = t_start;
                    while (t <= tf)
                    {
                        for (int i = 0; i < NJ; i++)
                        {
                            double pos, vel, acc;
                            calculateTrajectory(t, t_start, joint_coeffs[i], pos, vel, acc);
                            
                            POS_INIT(j * NJ + i) = pos;
                            VEL_INIT(j * NJ + i) = vel;
                            ACC_INIT(j * NJ + i) = acc;
                        }

                        std::vector<double> pos_des(NJ), vel_des(NJ), acc_des(NJ);
                        for(int k=0; k<NJ; k++) {
                             pos_des[k] = POS_INIT(j * NJ + k);
                             vel_des[k] = VEL_INIT(j * NJ + k);
                             acc_des[k] = ACC_INIT(j * NJ + k);
                        }

                        traj_msg.position = pos_des;
                        traj_msg.velocity = vel_des;
                        traj_msg.effort = acc_des;
                        
                        robot.set_q(Eigen::Map<Eigen::VectorXd>(pos_des.data(), NJ));
                        robot.set_dq(Eigen::Map<Eigen::VectorXd>(vel_des.data(), NJ));
                        robot.set_ddq(Eigen::Map<Eigen::VectorXd>(acc_des.data(), NJ));
                        
                        pub_cmd.publish(traj_msg);
                        publish_capsule_markers(robot, capsule_viz_pub, capsule_definitions, -1);

                        loop_rate_controller.sleep();
                        t = ros::Time::now().toSec();
                    }
                }

                // Invio pacchetto finale di stop (sicurezza)
                std::vector<double> final_pos(NJ), final_vel(NJ, 0.0), final_acc(NJ, 0.0);
                for(int i=0; i<NJ; i++) {
                    final_pos[i] = POS[(campioni - 1) * NJ + i];
                }
                traj_msg.position = final_pos;
                traj_msg.velocity = final_vel;
                traj_msg.effort = final_acc;
                pub_cmd.publish(traj_msg);

                std::cout << "Esecuzione completata." << std::endl;
            }
		}
	}
	return 0;
}