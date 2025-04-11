#include "utils/thunder_optimization.h"
#include <math.h>
#include <nlopt.hpp>
#include <stdio.h>

const std::string conf_file = "../robots/franka_conf.yaml";

std::vector<double> calculateCoefficients(double q0, double qf, double v0, double vf, double a0, double af, double t0, double tf)
{
    double T = tf - t0;

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
    double dt = t - t0;

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
void constraints(unsigned m, double *result, unsigned n, const double *x, double *grad, void *data)
{
    OptimizationData *optData = static_cast<OptimizationData *>(data);
    int NJ = optData->robot.get_numJoints();
    int campioni = optData->campioni;
    int size_q = optData->size_q;
    int z = 0;

    // std::cout << "Sono dentro i constarints:\n"<< optData->qf << std::endl;

    // Vincoli sulle condizioni iniziali
    for (int k = 0; k < NJ; k++)
    {
        result[z++] = x[k] - optData->q0[k];              // Posizione iniziale
        result[z++] = x[size_q + k] - optData->v0[k];     // Velocità iniziale
        result[z++] = x[2 * size_q + k] - optData->a0[k]; // Accelerazione iniziale
    }

    // Vincoli sulle condizioni finali
    for (int k = 0; k < NJ; k++)
    {
        result[z++] = x[(campioni - 1) * NJ + k] - optData->qf[k];
        result[z++] = x[size_q + (campioni - 1) * NJ + k] - optData->vf[k];
        result[z++] = x[2 * size_q + (campioni - 1) * NJ + k] - optData->af[k];
    }
    // for (int p=0; p <z ; p++)
    // {
    //     std::cout <<"result  "<< p << " - " << result[p]<< std::endl;
    // }
    // // Calcola tutti i vincoli
    // result[0] = x[7] - optData->qf[0];
    // result[1] = x[8] - optData->qf[1];
    // result[2] = x[9] - optData->qf[2];
    // result[3] = x[10] - optData->qf[3];
    // result[4] = x[11] - optData->qf[4];
    // result[5] = x[12] - optData->qf[5];
    // result[6] = x[13] - optData->qf[6];
    // result[7] = x[21] - optData->vf[0];
    // result[8] = x[22] - optData->vf[1];
    // result[9] = x[23] - optData->vf[2];
    // result[10] = x[24] - optData->vf[3];
    // result[11] = x[25] - optData->vf[4];
    // result[12] = x[26] - optData->vf[5];
    // result[13] = x[27] - optData->vf[6];
    // result[14] = x[35] - optData->af[0];
    // result[15] = x[36] - optData->af[1];
    // result[16] = x[37] - optData->af[2];
    // result[17] = x[38] - optData->af[3];
    // result[18] = x[39] - optData->af[4];
    // result[19] = x[40] - optData->af[5];
    // result[20] = x[41] - optData->af[6];

    if (grad != nullptr)
    {
        int z_vincolo = 0;
        std::fill_n(grad, m * n, 0.0);

        // Vincoli sulla condizione iniziale
        for (int k = 0; k < NJ; k++)
        {
            int index_q = k;
            grad[z_vincolo * n + index_q] = 1.0;
            z_vincolo++;
        }
        for (int k = 0; k < NJ; k++)
        {
            int index_dq = size_q + k;
            grad[z_vincolo * n + index_dq] = 1.0;
            z_vincolo++;
        }
        for (int k = 0; k < NJ; k++)
        {
            int index_ddq = 2 * size_q + k;
            grad[z_vincolo * n + index_ddq] = 1.0;
            z_vincolo++;
        }

        // Vincoli sulla condizione finale
        for (int k = 0; k < NJ; k++)
        {
            int index_q = (campioni - 1) * NJ + k;
            grad[z_vincolo * n + index_q] = 1.0;
            z_vincolo++;
        }
        for (int k = 0; k < NJ; k++)
        {
            int index_dq = size_q + (campioni - 1) * NJ + k;
            grad[z_vincolo * n + index_dq] = 1.0;
            z_vincolo++;
        }
        for (int k = 0; k < NJ; k++)
        {
            int index_ddq = 2 * size_q + (campioni - 1) * NJ + k;
            grad[z_vincolo * n + index_ddq] = 1.0;
            z_vincolo++;
        }
    }
}