
#include <math.h>
#include <nlopt.hpp>
#include <stdio.h>
#include <eigen3/Eigen/Dense>
#include "utils/thunder_franka.h"
#include "utils/thunder_optimization.h"
const std::string conf_file = "../config/franka_conf.yaml";
int main()
{
	thunder_franka robot;
	robot.load_conf(conf_file);
	const int NJ = robot.get_numJoints();
	// Parametri iniziali e finali per ciascun giunto

	Eigen::Vector<double, 7> q0(-1.25962, -0.663669, -0.692637, -2.17138, -0.264125, 1.50759, 0.0630972); // Posizioni iniziali
	Eigen::Vector<double, 7> qf(0.0, 0.0, 0.0, -0.2, 2.0, 1.0, 0.0);														   // Posizioni finali
	Eigen::Vector<double, 7> v0{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};															   // Velocità iniziali
	Eigen::Vector<double, 7> vf{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};															   // Velocità finali
	Eigen::Vector<double, 7> a0{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};															   // Accelerazioni iniziali
	Eigen::Vector<double, 7> af{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};															   // Accelerazioni finali

	double t0 = 0.0; // Tempo iniziale
	double tf = 10;	 // Tempo finale
	int time_step = 0;
	double dt = 0.1;				   // Passo di integrazione
	int intervalli = ((tf - t0) / dt); // Passo di tempo
	int campioni = (intervalli + 1);   // intervalli
	int size_q = NJ * campioni;		   // Dimensione di q
	// std::cout << "Campioni:" << campioni << std::endl;

	// Calcolo dei coefficienti per ciascun giunto
	std::vector<std::vector<double>> joint_coeffs(NJ);
	for (int i = 0; i < NJ; ++i)
	{
		joint_coeffs[i] = calculateCoefficients(q0[i], qf[i], v0[i], vf[i], a0[i], af[i], t0, tf);
	}

	// Generazione della traiettoria

	Eigen::MatrixXd M(NJ, NJ), C(NJ, NJ), G(NJ, 1);
	Eigen::MatrixXd POS(NJ, campioni), VEL(NJ, campioni), ACC(NJ, campioni), tau_minJerk(NJ, campioni);
	Eigen::VectorXd qa(NJ), dqa(NJ), ddqa(NJ), tauvec(NJ), tauveccol(NJ * campioni), POS_INIT(NJ * campioni), VEL_INIT(NJ * campioni), ACC_INIT(NJ * campioni);

	for (double t = t0; t <= tf; t += dt)
	{
		// std::cout << "Time: " << t << std::endl;
		for (int i = 0; i < NJ; i++)
		{
			double pos, vel, acc, tau;
			calculateTrajectory(t, t0, joint_coeffs[i], pos, vel, acc);
			// std::cout << "  Joint " << i << ": Pos = " << pos << ", Vel = " << vel << ", Acc = " << acc << std::endl;

			// qa(i) = pos;
			// dqa(i) = vel;
			// ddqa(i) = acc,
			POS(i, time_step) = pos;
			VEL(i, time_step) = vel;
			ACC(i, time_step) = acc;
			POS_INIT(time_step * NJ + i) = pos;
			VEL_INIT(time_step * NJ + i) = vel;
			ACC_INIT(time_step * NJ + i) = acc;
		}
		time_step++;
		std::cout << std::endl;
	}
	//  std::cout << "POS:" << POS_INIT << std::endl;
	//  std::cout << "VEL:" << VEL_INIT << std::endl;
	//  std::cout << "ACC:" << ACC_INIT << std::endl;
	time_step = 0;

	// Apertura del file CSV in modalità scrittura
	std::ofstream outfile1("risultato_min_jerk.csv");
	if (!outfile1)
	{
		std::cerr << "Errore nell'apertura del file per la scrittura." << std::endl;
		return 1;
	}
	// Scrittura dei dati: per ogni indice scrive una riga con q, dq e ddq separati da virgola
	for (size_t i = 0; i < size_q; ++i)
	{
		outfile1 << POS_INIT[i] << "," << VEL_INIT[i] << "," << ACC_INIT[i] << "\n";
	}

	outfile1.close();
	std::cout << "File CSV 'risultato_min_jerk.csv' scritto correttamente." << std::endl;

	for (int i = 0; i < campioni; ++i)
	{
		//   dqr = dqr.setOnes();
		//      ddqr = ddqr.setOnes();
		robot.set_q(POS.col(i));
		robot.set_dq(VEL.col(i));
		robot.set_ddq(ACC.col(i));
		// robot.setArguments(POS.col(i), VEL.col(i), dqr, ddqr);
		M = robot.get_M();
		C = robot.get_C();
		G = robot.get_G();
		// std::cout << "M:" << M << std::endl;
		// std::cout << "Pos col:" << i << POS.col(i) << std::endl;

		tauvec = M * ACC.col(i) + C * VEL.col(i) + G;
		// std::cout << "TAU:" << tauvec << std::endl;
		for (int j = 0; j < NJ; j++)
		{
			tau_minJerk(j, time_step) = tauvec(j);
			tauveccol(time_step * NJ + j) = tauvec(j);
		}

		time_step++;
	}

	// std::cout << "dim_tau_minJerk:" << tau_minJerk.size() << std::endl;
	// std::cout << "tau_minJerk:" << tau_minJerk << std::endl;
	//  std::cout << "tauveccol" << tauveccol<< std::endl;
	double costo_partenza = 0;
	int dim = tauveccol.size();
	std::cout << "dim:" << dim << std::endl;
	for (int i = 0; i < tauveccol.size(); i++)
	{
		costo_partenza += tauveccol[i] * tauveccol[i];
	};
	// std::cout<<"tauveccol:"<<tauveccol<<std::endl;
	std::cout << "costo_partenza:" << costo_partenza << std::endl;
	// Eigen::Vector<double, 7> qf1(0.0, 0.0, 0.0, 0.0, 7.0, 1.0, 0.0); // Posizioni finali

	OptimizationData optData;
	optData.robot; // Inizializza l'oggetto robot
	optData.q0 = q0;
	optData.v0 = v0;
	optData.a0 = a0;
	optData.qf = qf;
	optData.vf = vf;
	optData.af = af;
	optData.size_q = size_q;
	optData.dt = dt; // Passo temporale
	optData.campioni = campioni;

	// define the optimization problem
	nlopt::opt opt(nlopt::LD_MMA, 3 * NJ * campioni);
	opt.set_min_objective(objective, &optData);
	opt.set_xtol_rel(1e-4);

	std::vector<double> ub(3 * NJ * campioni), lb(3 * NJ * campioni), ubq(NJ), lbq(NJ), ubdq(NJ), lbdq(NJ), ubddq(NJ), lbddq(NJ);
	lbq = {-2.8973, -1.7628, -2.8973, -3.0718, -2.8973, -0.0175, -2.8973};
	ubq = {2.8973, 1.7628, 2.8973, -0.0698, 2.8973, 3.7525, 2.8973};
	lbdq = {-2.175, -2.175, -2.175, -2.175, -2.61, -2.61, -2.61};
	ubdq = {2.175, 2.175, 2.175, 2.175, 2.61, 2.61, 2.61};
	lbddq = {-15, -7.5, -10, -12.5, -15, -20, -20};
	ubddq = {15, 7.5, 10, 12.5, 15, 20, 20};
	
	float tol=1e-4;
//Vincoli sulle condizioni iniziali 
	for (int j = 0; j < 1; j++)
	{
		for (int i = 0; i < NJ; i++)
		{
			lb[j * NJ + i] = q0[i]-tol;
			ub[j * NJ + i] = q0[i]+tol;
			lb[j * NJ + size_q + i] = v0[i]-tol;
			ub[j * NJ + size_q + i] = v0[i]+tol;
			lb[j * NJ + 2 * size_q + i] = a0[i]-tol;
			ub[j * NJ + 2 * size_q + i] = a0[i]+tol;
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

//Vincoli sulle condizioni finali

	for (int j = campioni-1; j < campioni; j++)
	{
		for (int i = 0; i < NJ; i++)
		{
			lb[j * NJ + i] = qf[i]-tol;
			std::cout << "lb " << j * NJ + i << ":" << lb[j * NJ + i] << std::endl;
			ub[j * NJ + i] = qf[i]+tol;
			std::cout << "ub " << j * NJ + i << ":" << ub[j * NJ + i] << std::endl;
			lb[j * NJ + size_q + i] = vf[i]-tol;
			std::cout << "lb " << j * NJ + size_q + i << ":" << lb[j * NJ + size_q + i] << std::endl;
			ub[j * NJ + size_q + i] = vf[i]+tol;
			std::cout << "ub " << j * NJ + size_q + i << ":" << ub[j * NJ + size_q + i] << std::endl;
			lb[j * NJ + 2 * size_q + i] = af[i]-tol;
			std::cout << "lb " << j * NJ + 2*size_q + i << ":" << lb[j * NJ + 2*size_q + i] << std::endl;
			ub[j * NJ + 2 * size_q + i] = af[i]+tol;
			std::cout << "ub " << j * NJ + 2*size_q + i << ":" << ub[j * NJ + 2*size_q + i] << std::endl;
		}
	}
	opt.set_upper_bounds(ub);
	opt.set_lower_bounds(lb);

	std::vector<double> tol_constraint(NJ * 3 * 2);
	for (int i = 0; i < tol_constraint.size(); i++)
	{
		tol_constraint[i] = 1e-4;
	}

	opt.add_inequality_mconstraint(constraints, &optData, tol_constraint);

	// opt.set_maxtime(100);

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
	printf(",\n");
	printf("q_opt:\n");
	for (int i = 0; i < POS_INIT.size(); i++)
	{
		printf("%g ", vettore[i]);
		printf(",\n");
	}
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

	std::cout << "The result is " << result << std::endl;
	std::cout << "Minimal function value " << minf << std::endl;
	// Apertura del file CSV in modalità scrittura
	std::ofstream outfile("risultato_ottimo.csv");
	if (!outfile)
	{
		std::cerr << "Errore nell'apertura del file per la scrittura." << std::endl;
		return 1;
	}
	// Scrittura dei dati: per ogni indice scrive una riga con q, dq e ddq separati da virgola
	for (size_t i = 0; i < size_q; ++i)
	{
		outfile << vettore[i] << "," << vettore[size_q + i] << "," << vettore[2 * size_q + i] << "\n";
	}

	outfile.close();
	std::cout << "File CSV 'risultato_ottimo.csv' scritto correttamente." << std::endl;

	return 0;
}
