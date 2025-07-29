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

struct ObstacleConstraintIneq
{
    int k; // istante di tempo in cui valutare il vincolo
    int NJ;
    double r_s;
    double d_safe;
    Eigen::Vector3d p_obs;
    thunder_franka robot;
    Eigen::VectorXd q0;
    Eigen::VectorXd dq0;
    double dt;

    ObstacleConstraintIneq(int k_, int NJ_, double r_s_, double d_safe_,
                           const Eigen::Vector3d &p_obs_, thunder_franka robot_, Eigen::VectorXd q0_, Eigen::VectorXd dq0_, double dt_, Eigen::VectorXd qf_)
        : k(k_), NJ(NJ_), r_s(r_s_), d_safe(d_safe_), p_obs(p_obs_), robot(robot_), q0(q0_), dq0(dq0_), dt(dt_) {}
};

struct ConsistencyConstraintIneq
{
    int k;
    int NJ;
    int size_q;
    double dt;
    int type;     // 0 = posizione, 1 = velocità 2 = accelerazione
    int sign;     // +1 o -1 per la forma della disuguaglianza
    int campioni; // Numero di campioni per la traiettoria
    Eigen::VectorXd q0;
    Eigen::VectorXd v0;
    Eigen::VectorXd qf;
    int i; // indice del giunto su cui applicare il vincolo
};

struct LinkDistanceResult
{
    std::vector<double> distances;
    int closest_link_index;
    double min_distance;
};

struct JointLimitConstraint
{
    int NJ; // numero giunti
    int k;  // numero step di integrazione
    int i;  // indice giunto
    double dt;
    Eigen::VectorXd q0;
    Eigen::VectorXd v0;
    double limit;  // upper o lower limit
    bool is_upper; // true se è un upper bound, false se è un lower bound
};

// Funzione obiettivo per l'ottimizzazione
double objective(const std::vector<double> &x, std::vector<double> &grad, void *data);

// Funzione per i vincoli
double consistency_ineq(unsigned n, const double *x, double *grad, void *data);

// Funzione per evitare ostacoli sferici
double avoid_sphere(const std::vector<double> &x, std::vector<double> &grad, void *data);

// Funzione per evitare ostacoli sferici con gradiente
double avoid_sphere_with_gradient(const std::vector<double> &x, std::vector<double> &grad, void *data);

// Funzione per calcolare le distanze dei link da un punto (ostacolo)
LinkDistanceResult compute_link_distances_to_point(
    const Eigen::VectorXd &q,
    const Eigen::VectorXd &dq,
    const Eigen::VectorXd &ddq,
    const Eigen::Vector3d &p_obs,
    thunder_franka &robot);

// Vincolo per i limiti di posizione dei giunti
double joint_position_limit(unsigned n, const double *x, double *grad, void *data);

// Vincolo per i limiti di velocità dei giunti
double joint_velocity_limit(unsigned n, const double *x, double *grad, void *data);

#endif // THUNDER_OPTIMIZATION_H
