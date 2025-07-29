// Prima gli header standard C++
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// Poi gli header di sistema/ROS
#include "ros/ros.h"

// Poi le librerie di terze parti
#include <eigen3/Eigen/Dense>
#include <nlopt.hpp>

// Infine i tuoi header locali
#include "utils/thunder_optimization.h"

const std::string conf_file = "../robots/franka_conf.yaml";

std::vector<double> calculateCoefficients(double q0, double qf, double v0, double vf, double a0, double af, double t0, double tf)
{
    double T = tf;

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
    double dt = t;

    pos = coeffs[0] + coeffs[1] * dt + coeffs[2] * pow(dt, 2) + coeffs[3] * pow(dt, 3) + coeffs[4] * pow(dt, 4) + coeffs[5] * pow(dt, 5);
    vel = coeffs[1] + 2 * coeffs[2] * dt + 3 * coeffs[3] * pow(dt, 2) + 4 * coeffs[4] * pow(dt, 3) + 5 * coeffs[5] * pow(dt, 4);
    acc = 2 * coeffs[2] + 6 * coeffs[3] * dt + 12 * coeffs[4] * pow(dt, 2) + 20 * coeffs[5] * pow(dt, 3);
}

double objective(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    OptimizationData *optData = static_cast<OptimizationData *>(data);
    optData->robot.load_conf(conf_file);

    int NJ = optData->robot.get_numJoints();
    int campioni = optData->campioni;
    double dt = optData->dt;

    Eigen::VectorXd q = optData->q0;  // Inizializza q con q0
    Eigen::VectorXd dq = optData->v0; // Inizializza dq con v0
    Eigen::VectorXd ddq(NJ), tau_dyn(NJ);
    Eigen::VectorXd tau_tot(NJ * campioni); // salva tutte le tau

    double cost = 0.0;

    // Se gradiente richiesto, inizializza
    if (!grad.empty())
    {
        grad.assign(NJ * campioni, 0.0);
    }

    for (int k = 0; k < campioni; k++)
    {
        // estrai le ddq per il campione k
        for (int i = 0; i < NJ; i++)
        {
            ddq[i] = x[k * NJ + i];
            dq[i] += ddq[i] * dt;
            q[i] += dq[i] * dt;
        }

        // std::cout << "Campione: " << k << " q: " << q.transpose() << " dq: " << dq.transpose() << " ddq: " << ddq.transpose() << std::endl;

        // aggiorna stato robot
        optData->robot.set_q(q);
        optData->robot.set_dq(dq);
        optData->robot.set_ddq(ddq);

        Eigen::Matrix<double, 7, 7> M = optData->robot.get_M();
        Eigen::Matrix<double, 7, 7> C = optData->robot.get_C();
        Eigen::Matrix<double, 7, 1> G = optData->robot.get_G();

        tau_dyn = M * ddq + C * dq + G;

        // salva tau per cost finale
        for (int i = 0; i < NJ; i++)
        {
            tau_tot[k * NJ + i] = tau_dyn[i];
            cost += tau_dyn[i] * tau_dyn[i]; // somma al costo
        }

        // calcola gradiente rispetto a ddq_k
        if (!grad.empty())
        {
            Eigen::VectorXd grad_ddq = 2.0 * M.transpose() * tau_dyn;
            for (int i = 0; i < NJ; i++)
            {
                grad[k * NJ + i] = grad_ddq[i];
            }
        }
    }

    return cost;
}

// Vincolo per i limiti di posizione dei giunti

double joint_position_limit(unsigned n, const double *x, double *grad, void *data)
{
    JointLimitConstraint *c = reinterpret_cast<JointLimitConstraint *>(data);
    double dq = c->v0[c->i];
    double q = c->q0[c->i];

      // Se gradiente richiesto, inizializza
    if (grad)
    {
        std::fill(grad, grad + n, 0.0);
    }


    for (int j = 0; j < c->k; j++)
    {
        double ddq = x[j * c->NJ + c->i];
        dq += ddq * c->dt;
        q += dq * c->dt + 0.5 * ddq * c->dt * c->dt;
    }

    double val = c->is_upper ? (q - c->limit) : (c->limit - q);

    for (int j = 0; j < c->k; j++)
    {
        int idx = j * c->NJ + c->i;
        double coeff = (c->k - j) * c->dt * c->dt + 0.5 * c->dt * c->dt;
        grad[idx] = c->is_upper ? coeff : -coeff;
    }

    std::cout << "Valore vincolo posizione giunto " << c->i << ": " << val << std::endl;

    return val;
}

// Vincolo per i limiti di velocità dei giunti
double joint_velocity_limit(unsigned n, const double *x, double *grad, void *data)
{
    JointLimitConstraint *c = reinterpret_cast<JointLimitConstraint *>(data);
    double dq = c->v0[c->i];

    // Se gradiente richiesto, inizializza
    if (grad)
    {
        std::fill(grad, grad + n, 0.0);
    }

    for (int j = 0; j < c->k; j++)
    {
        double ddq = x[j * c->NJ + c->i];
        dq += ddq * c->dt;
    }

    double val = c->is_upper ? (dq - c->limit) : (c->limit - dq);

    for (int j = 0; j < c->k; j++)
    {
        int idx = j * c->NJ + c->i;
        grad[idx] = c->is_upper ? c->dt : -c->dt;
    }

    std::cout << "Valore vincolo velocità giunto " << c->i << ": " << val << std::endl;

    return val;
}

// // Funzione per i vincoli
double consistency_ineq(unsigned n, const double *x, double *grad, void *data)
{
    ConsistencyConstraintIneq *c = reinterpret_cast<ConsistencyConstraintIneq *>(data);

    double dq = c->v0[c->i];
    double q = c->q0[c->i];

    for (int j = 0; j < c->k; ++j)
    {
        double ddq = x[j * c->NJ + c->i];
        dq += ddq * c->dt;
        q += dq * c->dt + 0.5 * ddq * std::pow(c->dt, 2);
    }

    double err = q - c->qf[c->i];
    double val = err * err;

    if (grad)
    {
        std::fill(grad, grad + n, 0.0);
        for (int j = 0; j < c->k; ++j)
        {
            double coeff = (c->k - j) * c->dt * c->dt + 0.5 * c->dt * c->dt;
            int idx = j * c->NJ + c->i;
            grad[idx] = 2.0 * err * coeff;
        }
    }

    return val;
}

double avoid_sphere_with_gradient(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    ObstacleConstraintIneq *c = reinterpret_cast<ObstacleConstraintIneq *>(data);

    int k = c->k;
    int NJ = c->NJ;
    double r_s = c->r_s;
    double d_safe = c->d_safe;
    Eigen::Vector3d p_obs = c->p_obs;
    thunder_franka robot;
    // Ricostruzione di q_k e dq_k tramite integrazione (Eulero esplicito)
    Eigen::VectorXd q_k = c->q0;
    Eigen::VectorXd dq_k = c->dq0;
    double dt = c->dt;

    for (int j = 0; j < k; j++)
    {
        Eigen::VectorXd ddq_k(NJ);
        for (int i = 0; i < NJ; i++)
        {
            ddq_k[i] = x[j * NJ + i];
            dq_k[i] += ddq_k[i] * dt;                          // Aggiorna dq_k con la derivata
            q_k[i] += dq_k[i] * dt + 0.5 * ddq_k[i] * dt * dt; // Aggiorna q_k
        }
    }

    // Imposta configurazione del robot
    c->robot.set_q(q_k);

    // Calcola le posizioni dei link
    std::vector<Eigen::Vector3d> p_links = {
        c->robot.get_T_0_1().block<3, 1>(0, 3),
        c->robot.get_T_0_2().block<3, 1>(0, 3),
        c->robot.get_T_0_3().block<3, 1>(0, 3),
        c->robot.get_T_0_4().block<3, 1>(0, 3),
        c->robot.get_T_0_5().block<3, 1>(0, 3),
        c->robot.get_T_0_6().block<3, 1>(0, 3),
        c->robot.get_T_0_7().block<3, 1>(0, 3),
        c->robot.get_T_0_ee().block<3, 1>(0, 3) // End Effector
    };

    // Trova il link più vicino all'ostacolo
    double min_distance = 1e6;
    int closest_link = -1;
    Eigen::Vector3d closest_point;

    for (int i = 0; i < p_links.size(); ++i)
    {
        double distance = (p_links[i] - p_obs).norm();
        if (distance < min_distance)
        {
            min_distance = distance;
            closest_link = i;
            closest_point = p_links[i];
        }
    }

    std::cout << "Link più vicino: " << closest_link + 1 << " con distanza: " << min_distance << std::endl;

    // Calcola il valore del vincolo: deve essere <= 0
    double constraint_value = (r_s + d_safe) - min_distance;
    if (constraint_value < 0)
    {
        std::cout << "Vincolo soddisfatto: distanza sufficiente." << std::endl;
    }
    else
    {
        std::cout << "Vincolo non soddisfatto: distanza insufficiente." << std::endl;
    }

    std::cout << "Valore vincolo: " << constraint_value << std::endl;

    // Direzione di derivazione
    Eigen::Vector3d direction = (closest_point - p_obs).normalized();

    // Jacobiano del punto più vicino
    Eigen::MatrixXd J_closest(3, NJ);
    switch (closest_link)
    {
    case 0:
        J_closest = c->robot.get_J_1();
        break;
    case 1:
        J_closest = c->robot.get_J_2();
        break;
    case 2:
        J_closest = c->robot.get_J_3();
        break;
    case 3:
        J_closest = c->robot.get_J_4();
        break;
    case 4:
        J_closest = c->robot.get_J_5();
        break;
    case 5:
        J_closest = c->robot.get_J_6();
        break;
    case 6:
        J_closest = c->robot.get_J_7();
        break;
    case 7:
        J_closest = c->robot.get_J_ee();
        break;
    default:
        std::cerr << "Link non valido!" << std::endl;
        return constraint_value;
    }

    // Calcola derivata della distanza rispetto a q_k
    Eigen::RowVectorXd ddist_dq = direction.transpose() * J_closest;

    // Calcola derivata di q_k rispetto a tutte le ddq_j (j=0..k-1)
    // dq_k = dq_0 + ∑ ddq_j * dt  => d(dq_k)/d(ddq_j) = I * dt (per j < k)
    // q_k  = q_0  + ∑ dq_j * dt   => d(q_k)/d(ddq_j) = dt^2 * (k - j)

    // Reset gradiente
    if (!grad.empty())
        std::fill(grad.begin(), grad.end(), 0.0);

    for (int j = 0; j < k; ++j)
    {
        for (int i = 0; i < NJ; ++i)
        {
            int idx = j * NJ + i;
            double d_qk_i__d_ddqji = dt * dt * (k - j);
            grad[idx] = -ddist_dq[i] * d_qk_i__d_ddqji;
        }
    }

    return constraint_value;
}
