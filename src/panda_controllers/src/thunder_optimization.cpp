#include "utils/thunder_optimization.h"
#include <math.h>
#include <nlopt.hpp>
#include <stdio.h>
#include "ros/ros.h"

const std::string conf_file = "../robots/franka_conf.yaml";

std::vector<double> calculateCoefficients(double q0, double qf, double v0, double vf, double a0, double af, double t0, double tf)
{
    double T = tf ;

    double c0 = q0;
    double c1 = v0;
    double c2 = a0 / 2.0;
    double c3 = (20.0 * (qf - q0) - (8.0 * vf + 12.0 * v0) * T - (3.0 * a0 - af) * T * T) / (2.0 * pow(T, 3));
    double c4 = (30.0 * (q0 - qf) + (14.0 * vf + 16.0 * v0) * T + (3.0 * a0 - 2.0 * af) * T * T) / (2.0 * pow(T, 4));
    double c5 = (12.0 * (qf - q0) - 6.0 * (vf + v0) * T - (a0 - af) * T * T) / (2.0 * pow(T, 5));

    return {c0, c1, c2, c3, c4, c5};
}

void calculateTrajectory(double t, double t0, const std::vector<double> &coeffs, double &pos, double &vel, double &acc)
{
    double dt = t ;

    pos = coeffs[0] + coeffs[1] * dt + coeffs[2] * pow(dt, 2) + coeffs[3] * pow(dt, 3) + coeffs[4] * pow(dt, 4) + coeffs[5] * pow(dt, 5);
    vel = coeffs[1] + 2 * coeffs[2] * dt + 3 * coeffs[3] * pow(dt, 2) + 4 * coeffs[4] * pow(dt, 3) + 5 * coeffs[5] * pow(dt, 4);
    acc = 2 * coeffs[2] + 6 * coeffs[3] * dt + 12 * coeffs[4] * pow(dt, 2) + 20 * coeffs[5] * pow(dt, 3);
}

double objective(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    double cost = 0;
    OptimizationData *optData = static_cast<OptimizationData *>(data);
    optData->robot.load_conf(conf_file);
    int NJ = optData->robot.get_numJoints();
    int campioni = optData->campioni;
    int size_q = optData->size_q;

    Eigen::VectorXd tau_dyn(NJ), tau(NJ * campioni);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(NJ);
    Eigen::VectorXd dq = Eigen::VectorXd::Zero(NJ);
    Eigen::VectorXd ddq = Eigen::VectorXd::Zero(NJ);
    Eigen::VectorXd grad_q(NJ), grad_dq(NJ), grad_ddq(NJ), grade = Eigen::VectorXd::Zero(3 * NJ * campioni);

    for (int j = 0; j < optData->campioni; j++)
    {
        for (int i = 0; i < NJ; i++)
        {
            q[i] = x[j * NJ + i];
            dq[i] = x[j * NJ + size_q + i];
            ddq[i] = x[j * NJ + 2 * size_q + i];
        }

        optData->robot.set_q(q);
        optData->robot.set_dq(dq);
        optData->robot.set_ddq(ddq);

        Eigen::Matrix<double, 7, 7> M = optData->robot.get_M();
        Eigen::Matrix<double, 7, 7> C = optData->robot.get_C();
        Eigen::Matrix<double, 7, 1> G = optData->robot.get_G();

        tau_dyn = M * ddq + C * dq + G;

        for (int i = 0; i < NJ; i++)
        {
            tau[j * NJ + i] = tau_dyn[i];
        }

        grad_q = (2 * tau_dyn.transpose() * optData->robot.get_gradq());
        grad_dq = (2 * tau_dyn.transpose() * optData->robot.get_graddq());
        grad_ddq = (2 * tau_dyn.transpose() * optData->robot.get_gradddq());

        for (int k = 0; k < NJ; k++)
        {
            grade[j * NJ + k] = grad_q[k];
            grade[j * NJ + size_q + k] = grad_dq[k];
            grade[j * NJ + 2 * size_q + k] = grad_ddq[k];
        }

        if (!grad.empty())
        {
            for (int i = 0; i < 3 * NJ * campioni; i++)
            {
                grad[i] = grade[i];
            }
        }
    }

    for (int i = 0; i < NJ * campioni; i++)
    {
        cost += tau[i] * tau[i];
    }

    return cost;
}

// Funzione per i vincoli
double consistency_ineq(unsigned n, const double *x, double *grad, void *data) {
    ConsistencyConstraintIneq *c = reinterpret_cast<ConsistencyConstraintIneq *>(data);
    int k = c->k;
    int NJ = c->NJ;
    int size_q = c->size_q;
    double dt = c->dt;
    int type = c->type;
    int sgn = c->sign;

    double val = 0.0;
    for (int j = 0; j < NJ; ++j) {
        if (type == 0) {
            // Posizione: q_{k+1} - q_k - dq_k * dt
            std::cout <<"sono nel ciclo dei vincoli: "<< k << " type: "<< type << std::endl;
            int qk   = k * NJ + j;
            int qkp  = (k + 1) * NJ + j;
            int dqk  = size_q + k * NJ + j;
            val += sgn * (x[qkp] - x[qk] - x[dqk] * dt);
        } else if (type == 1) {
            // Velocità: dq_{k+1} - dq_k - ddq_k * dt
            std::cout <<"sono nel ciclo dei vincoli:  "<< k << " type: "<< type << std::endl;
            int dqk   = size_q + k * NJ + j;
            int dqkp  = size_q + (k + 1) * NJ + j;
            int ddqk  = 2 * size_q + k * NJ + j;
            val += sgn * (x[dqkp] - x[dqk] - x[ddqk] * dt);
        } else if (type == 2) {
            // Accelerazione: ddq_{k+1} - ddq_k
            std::cout <<"sono nel ciclo dei vincoli:  "<< k << " type: "<< type << std::endl;
            int ddqk   = 2 * size_q + k * NJ + j;
            int ddqkp  = 2 * size_q + (k + 1) * NJ + j;
            val += sgn * (x[ddqkp] - x[ddqk]);
        }
    }

    return val;
}
