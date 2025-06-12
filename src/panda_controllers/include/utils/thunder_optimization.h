#ifndef THUNDER_OPTIMIZATION_H
#define THUNDER_OPTIMIZATION_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include "thunder_franka.h"

// Dichiarazione della funzione per calcolare i coefficienti del polinomio di quinto grado
std::vector<double> calculateCoefficients(double q0, double qf, double v0, double vf, double a0, double af, double t0, double tf);

// Dichiarazione della funzione per calcolare posizione, velocità e accelerazione
void calculateTrajectory(double t, double t0, const std::vector<double> &coeffs, double &pos, double &vel, double &acc);

// Struttura per i dati di ottimizzazione
struct OptimizationData
{
    thunder_franka robot;
    Eigen::VectorXd q0, v0, a0, qf, vf, af;
    int size_q;
    double dt;
    int campioni;
};

struct ObstacleConstraintIneq {
    int k;
    int NJ;
    double r_s;
    double d_safe;
    Eigen::Vector3d p_obs;
    thunder_franka robot;

    ObstacleConstraintIneq(int k_, int NJ_, double r_s_, double d_safe_,
                           const Eigen::Vector3d& p_obs_, thunder_franka robot_)
        : k(k_), NJ(NJ_), r_s(r_s_), d_safe(d_safe_), p_obs(p_obs_), robot(robot_) {}
};


struct ConsistencyConstraintIneq
{
    int k;
    int NJ;
    int size_q;
    double dt;
    int type; // 0 = posizione, 1 = velocità 2 = accelerazione
    int sign; // +1 o -1 per la forma della disuguaglianza
};

// Funzione obiettivo per l'ottimizzazione
double objective(const std::vector<double> &x, std::vector<double> &grad, void *data);

// Funzione per i vincoli
double consistency_ineq(unsigned n, const double *x, double *grad, void *data);

// Funzione per evitare ostacoli sferici
double avoid_sphere(const std::vector<double> &x, std::vector<double> &grad, void *data);

#endif // THUNDER_OPTIMIZATION_H
