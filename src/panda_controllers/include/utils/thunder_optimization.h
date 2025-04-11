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

// Funzione obiettivo per l'ottimizzazione
double objective(const std::vector<double> &x, std::vector<double> &grad, void *data);

// Funzione per i vincoli
void constraints(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data);

#endif // THUNDER_OPTIMIZATION_H
