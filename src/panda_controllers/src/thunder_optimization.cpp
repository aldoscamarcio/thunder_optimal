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
double consistency_ineq(unsigned n, const double *x, double *grad, void *data)
{
    ConsistencyConstraintIneq *c = reinterpret_cast<ConsistencyConstraintIneq *>(data);
    int k = c->k;
    int NJ = c->NJ;
    int size_q = c->size_q;
    double dt = c->dt;
    int type = c->type;
    int sgn = c->sign;
    int campioni = c->campioni;

    for (int w = 0; w < 3 * NJ * campioni; w++)
    {
        grad[w] = 0.0; // Inizializzo il gradiente a zero
    }

    double val = 0.0;
    for (int j = 0; j < NJ; ++j)
    {
        if (type == 0)
        {
            // Posizione: q_{k+1} - q_k - dq_k * dt
            // std::cout <<"sono nel ciclo dei vincoli: "<< k << " type: "<< type << std::endl;
            int qk = k * NJ + j;
            int qkp = (k + 1) * NJ + j;
            int dqk = size_q + k * NJ + j;
            int ddqk = 2 * size_q + k * NJ + j;
          
            grad[qk] = -sgn * 1.0;
            grad[qkp] = sgn * 1.0;
            grad[dqk] = -sgn * dt;
            grad[ddqk] = -sgn * 0.5 * dt * dt;
            val += sgn * (x[qkp] - x[qk] - x[dqk] * dt -0.5*x[ddqk] * dt*dt);
        }
        else if (type == 1)
        {
            // Velocità: dq_{k+1} - dq_k - ddq_k * dt
            // std::cout <<"sono nel ciclo dei vincoli:  "<< k << " type: "<< type << std::endl;
            int dqk = size_q + k * NJ + j;
            int dqkp = size_q + (k + 1) * NJ + j;
            int ddqk = 2 * size_q + k * NJ + j;
            
            grad[dqk] = -sgn * 1.0;
            grad[dqkp] = sgn * 1.0;
            grad[ddqk] = -sgn * dt;
            // std::cout << "dqkp: " << dqkp << " dqk: " << dqk << " ddqk: " << ddqk << std::endl;
            val += sgn * (x[dqkp] - x[dqk] - x[ddqk] * dt);
        }
        else if (type == 2)
        {
            // Accelerazione: ddq_{k+1} - ddq_k
            // std::cout <<"sono nel ciclo dei vincoli:  "<< k << " type: "<< type << std::endl;
            int ddqk = 2 * size_q + k * NJ + j;
            int ddqkp = 2 * size_q + (k + 1) * NJ + j;
            grad[ddqk] = -sgn * 1.0;
            grad[ddqkp] = sgn * 1.0;
            // grad[ddqk] =  1.0;
            // std::cout << "ddqkp: " << ddqkp << " ddqk: " << ddqk << std::endl;
            val += sgn * (x[ddqkp] - x[ddqk]);
            // val += sgn * x[ddqk]-100;
        }
    }
    // std::cout << "valore vincolo: " << val << std::endl;
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

    // // Inizializza gradiente a zero
    // if (!grad.empty()) {
    //     std::fill(grad.begin(), grad.end(), 0.0);
    // }

    // Estrai q al tempo k
    int q_offset = k * NJ;
    Eigen::VectorXd q(NJ);
    for (int i = 0; i < NJ; i++) {
        q[i] = x[q_offset + i];
    }

    c->robot.set_q(q);

    // Calcola le posizioni dei link
    std::vector<Eigen::Vector3d> p_links = {
        c->robot.get_T_0_1().block<3,1>(0,3),
        c->robot.get_T_0_2().block<3,1>(0,3),
        c->robot.get_T_0_3().block<3,1>(0,3),
        c->robot.get_T_0_4().block<3,1>(0,3),
        c->robot.get_T_0_5().block<3,1>(0,3),
        c->robot.get_T_0_6().block<3,1>(0,3),
        c->robot.get_T_0_7().block<3,1>(0,3)
        //c->robot.get_T_0_ee().block<3,1>(0,3)   senza gripper
    };

    // Trova il link più vicino
    double min_distance = 1e6;
    int closest_link = -1;
    Eigen::Vector3d closest_point;

    for (int i = 0; i < p_links.size(); ++i) {
        double distance = (p_links[i] - p_obs).norm();
        if (distance < min_distance) {
            min_distance = distance;
            closest_link = i;
            closest_point = p_links[i];
        }
    }

    std::cout << "Link più vicino: " << closest_link + 1 << std::endl;
    std::cout << "Distanza minima: " << min_distance << std::endl;

    // Vincolo di disuguaglianza: deve essere <= 0
    double constraint_value = (r_s + d_safe) - min_distance;
    std::cout << "constraint value: " << constraint_value << std::endl;

        Eigen::Vector3d direction = (closest_point - p_obs).normalized();

        // Recupera lo Jacobiano giusto
        Eigen::MatrixXd J_closest(3, NJ);
        switch (closest_link) {
            case 0: J_closest = c->robot.get_J_1(); break;
            case 1: J_closest = c->robot.get_J_2(); break;
            case 2: J_closest = c->robot.get_J_3(); break;
            case 3: J_closest = c->robot.get_J_4(); break;
            case 4: J_closest = c->robot.get_J_5(); break;
            case 5: J_closest = c->robot.get_J_6(); break;
            case 6: J_closest = c->robot.get_J_7(); break;
            case 7: J_closest = c->robot.get_J_ee(); break;
            default:
                std::cerr << "Link non valido per Jacobiano!" << std::endl;
                return constraint_value;
        }

        // Calcolo derivata della distanza rispetto a q
        Eigen::VectorXd ddist_dq = direction.transpose() * J_closest;

        // Inserisco nel vettore grad
        for (int i = 0; i < NJ; i++) {
            grad[q_offset + i] = -ddist_dq[i];  // -ddist_dq perché constraint = (r_s + d_safe) - dist
        }


    return constraint_value;
}

//funzione per calcolare le distanze dei link da un punto (ostacolo)
LinkDistanceResult compute_link_distances_to_point(
    const Eigen::VectorXd& q,
    const Eigen::VectorXd& dq,
    const Eigen::VectorXd& ddq,
    const Eigen::Vector3d& p_obs,
    thunder_franka& robot
) {
    robot.set_q(q);
    robot.set_dq(dq);
    robot.set_ddq(ddq);

    std::vector<Eigen::Vector3d> p_links = {
        robot.get_T_0_1().block<3,1>(0,3),
        robot.get_T_0_2().block<3,1>(0,3),
        robot.get_T_0_3().block<3,1>(0,3),
        robot.get_T_0_4().block<3,1>(0,3),
        robot.get_T_0_5().block<3,1>(0,3),
        robot.get_T_0_6().block<3,1>(0,3),
        robot.get_T_0_7().block<3,1>(0,3)
        // robot.get_T_0_ee().block<3,1>(0,3)
    };

    std::vector<double> distances;
    distances.reserve(p_links.size());

    double min_distance = std::numeric_limits<double>::max();
    int closest_index = -1;

    for (size_t i = 0; i < p_links.size(); ++i) {
        double d = (p_links[i] - p_obs).norm();
        distances.push_back(d);
        if (d < min_distance) {
            min_distance = d;
            closest_index = static_cast<int>(i);
        }
    }

    return {distances, closest_index, min_distance};
}
