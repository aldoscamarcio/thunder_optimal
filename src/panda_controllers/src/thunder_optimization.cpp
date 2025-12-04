// Prima gli header standard C++
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

// Poi gli header di sistema/ROS
#include "ros/ros.h"
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <sstream> // ✅ Per std::stringstream
#include <iomanip>

// Poi le librerie di terze parti
#include <eigen3/Eigen/Dense>
#include <nlopt.hpp>

// Infine i tuoi header locali
#include "utils/thunder_optimization.h"
#include <geometry_msgs/Point.h>

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

    // std::cout << "Valore vincolo posizione giunto " << c->i << ": " << val << std::endl;

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

    // std::cout << "Valore vincolo velocità giunto " << c->i << ": " << val << std::endl;

    return val;
}

// Vincolo sulla posizione finale del giunto i
double final_position_constraint(unsigned n, const double *x, double *grad, void *data)
{
    ConsistencyConstraintIneq *c = reinterpret_cast<ConsistencyConstraintIneq *>(data);

    // ATTENZIONE! k = campioni = cost; quindi lui ogni volta che ottimizza integra fino alla qf e la confronta con la qf settata

    double dq = c->v0[c->i];
    double q = c->q0[c->i];

    for (int j = 0; j < c->k; ++j)
    {
        double ddq = x[j * c->NJ + c->i];
        dq += ddq * c->dt;
        q += dq * c->dt + 0.5 * ddq * std::pow(c->dt, 2);
    }
    // std::cout << "i: " << c->i << " q: " << q << std::endl;

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

// Vincolo "soft" sulla velocità finale del giunto i
double final_velocity_constraint(unsigned n, const double *x, double *grad, void *data)
{
    ConsistencyConstraintIneq *c = reinterpret_cast<ConsistencyConstraintIneq *>(data);

    double dq = c->v0[c->i]; // velocità iniziale

    // Ricostruzione della velocità finale tramite integrazione delle accelerazioni
    for (int j = 0; j < c->k; ++j)
    {
        double ddq = x[j * c->NJ + c->i];
        dq += ddq * c->dt;
    }
    // std::cout << "i: " << c->i << " dq: " << dq << std::endl;

    // Valore del vincolo quadratico (soft constraint)
    double err = dq - c->vf[c->i]; // differenza tra velocità finale e target
    double val = err * err;

    // Gradiente
    if (grad)
    {
        std::fill(grad, grad + n, 0.0);
        for (int j = 0; j < c->k; ++j)
        {
            int idx = j * c->NJ + c->i;
            grad[idx] = 2.0 * err * c->dt; // derivata del quadrato rispetto a ddq_j
        }
    }

    return val;
}

void publish_capsule_markers(
    thunder_franka &robot,                            // Robot con q impostato
    ros::Publisher &marker_pub,                       // Publisher (passato per riferimento)
    const std::vector<Capsule> &capsules_definitions, // Definizioni delle capsule
    int closest_capsule_index)                        // per colorare la più vicina
{
    // 1. Ottieni le pose (basate sullo stato 'q' già impostato nel robot)
    std::vector<Eigen::Matrix4d> link_poses = {
        robot.get_T_0_0(), // Indice 0
        robot.get_T_0_1(), // Indice 1
        robot.get_T_0_2(), // Indice 2
        robot.get_T_0_3(), // Indice 3
        robot.get_T_0_4(), // Indice 4
        robot.get_T_0_5(), // Indice 5
        robot.get_T_0_5(), // Indice 5
        robot.get_T_0_6(), // Indice 6
        robot.get_T_0_7()  // Indice 7
    };

    // 2. Crea l'array di marker
    visualization_msgs::MarkerArray marker_array;

    for (size_t i = 0; i < capsules_definitions.size(); ++i)
    {
        const auto &cap = capsules_definitions[i];

        if (cap.link_index < 0 || cap.link_index >= link_poses.size())
            continue;

        const Eigen::Matrix4d &T_world_link = link_poses[cap.link_index];
        Eigen::Matrix4d T_world_capsule = T_world_link * cap.T_offset;

        Eigen::Vector3d center = T_world_capsule.block<3, 1>(0, 3);
        Eigen::Vector3d z_axis = T_world_capsule.block<3, 1>(0, 2);
        Eigen::Vector3d half_axis = z_axis * (cap.length / 2.0);
        Eigen::Vector3d a = center - half_axis;
        Eigen::Vector3d b = center + half_axis;

        // 3. Crea il Marker
        visualization_msgs::Marker marker;
        marker.header.frame_id = "panda_link0";
        marker.header.stamp = ros::Time::now();
        marker.ns = "collision_capsules";
        marker.id = static_cast<int>(i);
        marker.type = visualization_msgs::Marker::CYLINDER;
        marker.action = visualization_msgs::Marker::ADD;

        // 4. Posa (Centro + Orientamento)
        Eigen::Vector3d marker_center = (a + b) / 2.0;
        Eigen::Vector3d axis_vector = (b - a).normalized();
        Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), axis_vector);

        marker.pose.position.x = marker_center.x();
        marker.pose.position.y = marker_center.y();
        marker.pose.position.z = marker_center.z();
        marker.pose.orientation.x = q.x();
        marker.pose.orientation.y = q.y();
        marker.pose.orientation.z = q.z();
        marker.pose.orientation.w = q.w();

        // 5. Scala
        marker.scale.x = cap.radius * 2.0;
        marker.scale.y = cap.radius * 2.0;
        marker.scale.z = cap.length;

        // 6. Colore
        marker.color.r = 0.0f;
        marker.color.g = 1.0f;
        marker.color.b = 0.0f;
        marker.color.a = 0.4f;
        if (static_cast<int>(i) == closest_capsule_index)
        {
            marker.color.r = 1.0f;
            marker.color.g = 0.0f;
            marker.color.a = 0.8f;
        }
        marker.lifetime = ros::Duration(2.0);

        marker_array.markers.push_back(marker);
    }

    // 4. Pubblica
    marker_pub.publish(marker_array);
}

// Funzione per distanza punto-capsula
double point_to_capsule_distance(const Eigen::Vector3d &p, const Eigen::Vector3d &a, const Eigen::Vector3d &b, double radius)
{
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ap = p - a;
    double t = ap.dot(ab) / ab.squaredNorm(); // semplificazione
    t = std::min(std::max(t, 0.0), 1.0);
    Eigen::Vector3d closest = a + t * ab;
    // std::cout << "Distance: " << (p - closest).norm() - radius << std::endl;
    return (p - closest).norm() - radius;
}

// // Funzione di vincolo con gradiente
// double avoid_sphere_with_gradient(const std::vector<double> &x, std::vector<double> &grad, void *data)
// {
//     ObstacleConstraintIneq *c = reinterpret_cast<ObstacleConstraintIneq *>(data);

//     int k = c->k;
//     int NJ = c->NJ;
//     double r_s = c->r_s;
//     double d_safe = c->d_safe;
//     Eigen::Vector3d p_obs = c->p_obs;
//     double dt = c->dt;

//     Eigen::VectorXd q_k = c->q0;
//     Eigen::VectorXd dq_k = c->dq0;
//     for (int j = 0; j < k; j++)
//     {
//         for (int i = 0; i < NJ; i++)
//         {
//             double ddq = x[j * NJ + i];
//             dq_k[i] += ddq * dt;
//             q_k[i] += dq_k[i] * dt;
//         }
//     }

//     // Aggiorna configurazione robot
//     c->robot.set_q(q_k);

//     // Costruisce pose e Jacobiani di tutti i link
//     std::vector<Eigen::Matrix4d> link_poses = {
//         c->robot.get_T_0_0(),
//         c->robot.get_T_0_1(),
//         c->robot.get_T_0_2(),
//         c->robot.get_T_0_3(),
//         c->robot.get_T_0_4(),
//         c->robot.get_T_0_5(),
//         c->robot.get_T_0_5(),
//         c->robot.get_T_0_6(),
//         c->robot.get_T_0_7()};

//     std::vector<Eigen::MatrixXd> J_links = {
//         Eigen::MatrixXd::Zero(6, 7),
//         c->robot.get_J_1(),
//         c->robot.get_J_2(),
//         c->robot.get_J_3(),
//         c->robot.get_J_4(),
//         c->robot.get_J_5(),
//         c->robot.get_J_5(),
//         c->robot.get_J_6(),
//         c->robot.get_J_7()};

//     // Trova la capsula più vicina all'ostacolo
//     double min_distance = 1e6;
//     Eigen::Vector3d closest_point;
//     Eigen::MatrixXd J_closest(3, NJ);
//     J_closest.setZero();

//     int closest_capsule_index = -1;
//     for (size_t idx = 0; idx < c->capsules.size(); idx++)
//     {
//         auto &cap = c->capsules[idx];
//         std::cout << "numero capsule"<< c->capsules.size()<< std::endl;
//         Eigen::Matrix4d T = link_poses[cap.link_index] * cap.T_offset;
//         Eigen::Vector3d a = T.block<3, 1>(0, 3);
//         Eigen::Vector3d b = a + T.block<3, 1>(0, 2) * cap.length;
//         std::cout << "Lunghezza Capsula " << cap.length << std::endl;
//         double dist = point_to_capsule_distance(p_obs, a, b, cap.radius);

//         if (dist < min_distance)
//         {
//             min_distance = dist;
//             Eigen::Vector3d ab = b - a;
//             double t = ((p_obs - a).dot(ab)) / ab.squaredNorm();
//             t = std::min(std::max(t, 0.0), 1.0);
//             closest_point = a + t * ab;
//             closest_capsule_index = idx; // Salva indice

//             // 1. Prendi Jacobiani del LINK
//             Eigen::MatrixXd J_link_trans = J_links[cap.link_index].block(0, 0, 3, NJ);
//             Eigen::MatrixXd J_link_rot = J_links[cap.link_index].block(3, 0, 3, NJ);

//             // 2. Posizione del frame del link
//             Eigen::Vector3d p_link = link_poses[cap.link_index].block<3, 1>(0, 3);

//             // 3. Vettori offset (in coordinate globali)
//             Eigen::Vector3d r_link_to_a = a - p_link;
//             Eigen::Vector3d r_link_to_b = b - p_link;

//             // 4. Calcola J_a e J_b
//             Eigen::MatrixXd J_a(3, NJ);
//             Eigen::MatrixXd J_b(3, NJ);
//             for (int col = 0; col < NJ; ++col)
//             {
//                 Eigen::Vector3d omega = J_link_rot.col(col);

//                 // omega per il prodotto vettoriale
//                 J_a.col(col) = J_link_trans.col(col) + omega.cross(r_link_to_a);
//                 J_b.col(col) = J_link_trans.col(col) + omega.cross(r_link_to_b);
//             }

//             // 5. Interpola per trovare J_closest
//             J_closest = (1.0 - t) * J_a + t * J_b;
//         }
//     }

//     if (closest_capsule_index == -1)
//     {
//         return 0.0; // O Houston, abbiamo un problema
//     }

//     std::cout << "Capsula più vicina: " << closest_capsule_index << std::endl;

//     // Vincolo = distanza minima - sicurezza
//     double constraint_value = (r_s + d_safe) - min_distance;
//     std::cout << "Distanza minima: " << constraint_value << std::endl;
//     if (constraint_value < 0)
//     {
//         std::cout << "Nessun rischio di collisione." << std::endl;
//     }

//     // Visualizzazione
//     if (c->marker_pub)
//     {
//         publish_capsule_markers(c->robot, c->marker_pub, c->capsules, closest_capsule_index);
//     }

//     // Gradiente
//     Eigen::Vector3d direction = (closest_point - p_obs).normalized();
//     Eigen::RowVectorXd ddist_dq = direction.transpose() * J_closest;

//     if (!grad.empty())
//         std::fill(grad.begin(), grad.end(), 0.0);

//     // Loop del gradiente
//     for (int j = 0; j < k; ++j)
//     {
//         for (int i = 0; i < NJ; ++i)
//         {
//             int idx = j * NJ + i;
//             double d_qk_i__d_ddqji = dt * dt * (k - j);
//             grad[idx] = -ddist_dq[i] * d_qk_i__d_ddqji;
//         }
//     }

//     return constraint_value;
// }
double avoid_obstacle_generic(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    ObstacleConstraintIneq *c = reinterpret_cast<ObstacleConstraintIneq *>(data);

    int k = c->k;
    int NJ = c->NJ;
    double dt = c->dt;

    // 1. Integrazione in avanti per ottenere q(k) dalle variabili di decisione (ddq)
    Eigen::VectorXd q_k = c->q0;
    Eigen::VectorXd dq_k = c->dq0;
    for (int j = 0; j < k; j++)
    {
        for (int i = 0; i < NJ; i++)
        {
            double ddq = x[j * NJ + i];
            dq_k[i] += ddq * dt;
            q_k[i] += dq_k[i] * dt;
        }
    }

    // 2. Aggiorna cinematica robot
    c->robot->set_q(q_k);

    // Recupera pose e Jacobiani (come nel tuo script originale)
    // Nota: Assumo che get_T... e get_J... siano metodi della tua classe Robot
    std::vector<Eigen::Matrix4d> link_poses = {
        c->robot->get_T_0_0(), c->robot->get_T_0_1(), c->robot->get_T_0_2(), c->robot->get_T_0_3(),
        c->robot->get_T_0_4(), c->robot->get_T_0_5(), c->robot->get_T_0_5(), /*link flange?*/
        c->robot->get_T_0_6(), c->robot->get_T_0_7()};

    std::vector<Eigen::MatrixXd> J_links = {
        Eigen::MatrixXd::Zero(6, 7), c->robot->get_J_1(), c->robot->get_J_2(), c->robot->get_J_3(),
        c->robot->get_J_4(), c->robot->get_J_5(), c->robot->get_J_5(),
        c->robot->get_J_6(), c->robot->get_J_7()};

    // 3. Trova la distanza minima tra TUTTE le capsule del robot e l'OSTACOLO
    double min_signed_dist = 1e6;
    Eigen::MatrixXd J_closest(3, NJ);
    J_closest.setZero();
    Eigen::Vector3d gradient_direction = Eigen::Vector3d::Zero();

    int closest_cap_idx = -1;

    for (size_t idx = 0; idx < c->capsules_definitions.size(); idx++)
    {
        auto &cap_def = c->capsules_definitions[idx];

        // Calcola A e B in World Frame
        Eigen::Matrix4d T = link_poses[cap_def.link_index] * cap_def.T_offset;
        Eigen::Vector3d A_world = T.block<3, 1>(0, 3);
        // Nota: Nel tuo URDF le capsule sembrano allineate lungo Z locale
        Eigen::Vector3d B_world = A_world + T.block<3, 1>(0, 2) * cap_def.length;

        // Costruisci la capsula per il framework DistanceFunctions
        CapsuleWorld cap_world;
        cap_world.A = A_world;
        cap_world.B = B_world;
        cap_world.radius = cap_def.radius;

        // --- CHIAMATA AL NUOVO FRAMEWORK ---
        CapsuleDistanceResult res;
        double dist = capsule_distance(cap_world, c->obstacle, &res);
        // -----------------------------------

        if (dist < min_signed_dist)
        {
            min_signed_dist = dist;
            closest_cap_idx = idx;

            // La normale 'res.normal' punta DA ostacolo A capsula robot.
            // Per allontanarci, dobbiamo muoverci lungo questa normale.
            gradient_direction = res.normal;

            // --- Calcolo Jacobiano nel punto di minima distanza ---
            // Usiamo 'res.t_capsule' che ci dice esattamente dove cade il punto sul segmento
            double t = res.t_capsule; // 0.0 = A, 1.0 = B

            int link_idx = cap_def.link_index;
            Eigen::MatrixXd J_trans = J_links[link_idx].block(0, 0, 3, NJ);
            Eigen::MatrixXd J_rot = J_links[link_idx].block(3, 0, 3, NJ);
            Eigen::Vector3d p_link_origin = link_poses[link_idx].block<3, 1>(0, 3);

            // Bracci di leva per A e B
            Eigen::Vector3d r_A = A_world - p_link_origin;
            Eigen::Vector3d r_B = B_world - p_link_origin;

            // Jacobiano geometrico interpolato
            Eigen::MatrixXd J_A(3, NJ), J_B(3, NJ);
            for (int col = 0; col < NJ; ++col)
            {
                Eigen::Vector3d omega = J_rot.col(col);
                J_A.col(col) = J_trans.col(col) + omega.cross(r_A);
                J_B.col(col) = J_trans.col(col) + omega.cross(r_B);
            }

            // Interpolazione lineare precisa grazie a t_capsule
            J_closest = (1.0 - t) * J_A + t * J_B;
        }
    }

    // if (closest_cap_idx != -1 && k > 90)
    // {
    //     std::cout << "[Step " << k << "] Closest Capsule ID: " << closest_cap_idx
    //               << " | Dist: " << min_signed_dist << std::endl;
    // }

    if (closest_cap_idx == -1)
        return -1.0; // Fallback

    // 4. Definizione del valore del vincolo
    // Vogliamo: dist > d_safe  =>  d_safe - dist < 0
    // min_signed_dist è la distanza effettiva tra le superfici (già sottratti i raggi)
    double constraint_value = c->d_safe - min_signed_dist;

    /*
    std::cout << "Min Dist: " << min_signed_dist
              << " | Cons: " << constraint_value
              << " | Obs Type: " << (int)c->obstacle.type << std::endl;
    */

    // Visualizzazione
    if (c->marker_pub)
    {
        publish_capsule_markers(*c->robot, *c->marker_pub, c->capsules_definitions, closest_cap_idx);
    }

    // 5. Calcolo del Gradiente Analitico per NLopt
    if (!grad.empty())
    {
        // Chain rule: d(Constraint)/dq = - d(dist)/dq
        // d(dist)/dq = n_transpose * J_point
        // Quindi: Grad_C = - (n_transpose * J)
        Eigen::RowVectorXd ddist_dq = gradient_direction.transpose() * J_closest;
        Eigen::RowVectorXd dConstraint_dq = -ddist_dq; // Il segno meno è perché C = safe - dist

        // Propagazione indietro nel tempo (da q a ddq)
        // d_qk / d_ddq = dt^2 * (k - j)
        for (int j = 0; j < k; ++j)
        {
            double scaling_factor = dt * dt * (k - j);
            for (int i = 0; i < NJ; ++i)
            {
                grad[j * NJ + i] = dConstraint_dq[i] * scaling_factor;
            }
        }

        // Per i passi j >= k il gradiente è 0 (causalità)
        std::fill(grad.begin() + k * NJ, grad.end(), 0.0);
    }

    return constraint_value;
}

double avoid_self_collision(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    SelfCollisionConstraint *c_self = reinterpret_cast<SelfCollisionConstraint *>(data);

    int k = c_self->k;
    int NJ = c_self->NJ;
    double dt = c_self->dt;

    // 1. Integrazione per ottenere q(k)
    Eigen::VectorXd q_k = c_self->q0;
    Eigen::VectorXd dq_k = c_self->dq0;

    for (int j = 0; j < k; j++)
    {
        for (int i = 0; i < NJ; i++)
        {
            double ddq = x[j * NJ + i];

            dq_k(i) += ddq * dt;                          // v(t+1) = v(t) + a*dt
            q_k(i) += dq_k(i) * dt + 0.5 * ddq * dt * dt; // q(t+1) = q(t) + v*dt + 0.5*a*dt^2
        }
    }

    // 2. Aggiorna cinematica
    c_self->robot->set_q(q_k);

    // 3. Recupera pose e Jacobiani per tutti i link
    std::vector<Eigen::Matrix4d> link_poses = {
        c_self->robot->get_T_0_0(), c_self->robot->get_T_0_0(), c_self->robot->get_T_0_1(), c_self->robot->get_T_0_2(),
        c_self->robot->get_T_0_3(), c_self->robot->get_T_0_4(), c_self->robot->get_T_0_5(), c_self->robot->get_T_0_5(),
        c_self->robot->get_T_0_6(), c_self->robot->get_T_0_7(), c_self->robot->get_T_0_7(), c_self->robot->get_T_0_8()};

    std::vector<Eigen::MatrixXd> J_links = {
        Eigen::MatrixXd::Zero(6, NJ), // Link 0 (base fissa)
        c_self->robot->get_J_1(), c_self->robot->get_J_2(), c_self->robot->get_J_3(),
        c_self->robot->get_J_4(), c_self->robot->get_J_5(), c_self->robot->get_J_6(),
        c_self->robot->get_J_7(), c_self->robot->get_J_ee()};

    // 4. Calcola posizioni world di tutte le capsule
    std::vector<CapsuleWorld> capsules_world;
    capsules_world.reserve(c_self->capsules_definitions.size());

    for (const auto &cap_def : c_self->capsules_definitions)
    {
        Eigen::Matrix4d T = link_poses[cap_def.link_index] * cap_def.T_offset;

        CapsuleWorld cap_world;
        cap_world.A = T.block<3, 1>(0, 3);
        cap_world.B = cap_world.A + T.block<3, 1>(0, 2) * cap_def.length;
        cap_world.radius = cap_def.radius;

        capsules_world.push_back(cap_world);
    }

    // 5. Trova distanza minima tra coppie di link configurate
    double min_signed_dist = 1e6;
    int closest_cap_i = -1;
    int closest_cap_j = -1;
    CapsuleDistanceResult best_result;

    // Itera su tutte le coppie di collision configurate
    for (const auto &pair : c_self->collision_pairs)
    {
        int link_i = pair.first;
        int link_j = pair.second;

        // Trova tutte le capsule appartenenti a link_i e link_j
        for (size_t idx_i = 0; idx_i < c_self->capsules_definitions.size(); idx_i++)
        {
            if (c_self->capsules_definitions[idx_i].link_index != link_i)
                continue;

            for (size_t idx_j = 0; idx_j < c_self->capsules_definitions.size(); idx_j++)
            {
                if (c_self->capsules_definitions[idx_j].link_index != link_j)
                    continue;

                // Calcola distanza tra le due capsule
                CapsuleDistanceResult res;
                double dist = dist_capsule_capsule(capsules_world[idx_i],
                                                   capsules_world[idx_j],
                                                   &res);

                if (dist < min_signed_dist)
                {
                    min_signed_dist = dist;
                    closest_cap_i = idx_i;
                    closest_cap_j = idx_j;
                    best_result = res;
                }
            }
        }
    }

    // Nessuna coppia trovata (non dovrebbe succedere)
    if (closest_cap_i == -1 || closest_cap_j == -1)
    {
        if (!grad.empty())
            std::fill(grad.begin(), grad.end(), 0.0);
        return -1.0;
    }

    // 6. Definizione vincolo: dist > d_safe => d_safe - dist < 0
    double constraint_value = c_self->d_safe - min_signed_dist;

    static int debug_calls = 0;
    debug_calls++;
    if (debug_calls % 1 == 0)
    {
        std::cout << "[k=" << c_self->k << ", call #" << debug_calls << "]"
                  << " dist=" << min_signed_dist
                  << ", d_safe=" << c_self->d_safe
                  << ", constraint=" << constraint_value
                  << " | Caps " << closest_cap_i << " (Link "
                  << c_self->capsules_definitions[closest_cap_i].link_index
                  << ") <-> " << closest_cap_j << " (Link "
                  << c_self->capsules_definitions[closest_cap_j].link_index << ")"
                  << std::endl;
    }

    // Debug output (throttled)
    if (k % 10 == 0 && constraint_value > -0.05)
    {
        ROS_DEBUG("Step %d: Self-collision dist=%.3f (caps %d<->%d, links %d<->%d)",
                  k, min_signed_dist, closest_cap_i, closest_cap_j,
                  c_self->capsules_definitions[closest_cap_i].link_index,
                  c_self->capsules_definitions[closest_cap_j].link_index);
    }

    // 7. Calcolo Gradiente Analitico
    if (!grad.empty())
    {
        std::fill(grad.begin(), grad.end(), 0.0);

        // Normale di collisione: punta da capsula_j verso capsula_i
        Eigen::Vector3d normal = best_result.normal;

        // Parametri per interpolazione sui segmenti
        double t_i = best_result.t_capsule;  // Parametro su capsula_i
        double t_j = best_result.t_obstacle; // Parametro su capsula_j

        // Recupera info delle due capsule
        auto &cap_i_def = c_self->capsules_definitions[closest_cap_i];
        auto &cap_j_def = c_self->capsules_definitions[closest_cap_j];

        int link_i = cap_i_def.link_index;
        int link_j = cap_j_def.link_index;

        // Jacobiano per capsula i
        Eigen::MatrixXd J_i = compute_capsule_jacobian(
            link_poses[link_i], cap_i_def, J_links[link_i], t_i, NJ);

        // Jacobiano per capsula j
        Eigen::MatrixXd J_j = compute_capsule_jacobian(
            link_poses[link_j], cap_j_def, J_links[link_j], t_j, NJ);

        // Gradiente della distanza rispetto a q:
        // d(dist)/dq = normal^T * (J_i - J_j)
        Eigen::RowVectorXd ddist_dq = normal.transpose() * (J_i - J_j);

        // Gradiente del vincolo: d(C)/dq = -d(dist)/dq
        Eigen::RowVectorXd dConstraint_dq = -ddist_dq;

        for (int j = 0; j < k; j++)
        {
            // Chain rule: dC/d(ddq_j) = dC/dq_k * dq_k/d(ddq_j)
            // dq_k/d(ddq_j) = sum_{t=j}^{k-1} dt^2 * (1 + (t-j))

            for (int i = 0; i < NJ; i++)
            {
                double dq_k_dddq = 0.0;

                // Integra contributo da timestep j fino a k
                for (int t = j; t < k; t++)
                {
                    dq_k_dddq += dt * dt * (0.5 + (t - j));
                }

                grad[j * NJ + i] = dConstraint_dq[i] * dq_k_dddq;
            }
        }
    }

    // 8. Visualizzazione (opzionale)
    if (c_self->marker_pub && (c_self->k == 25 || constraint_value > -0.02))
    {
        publish_self_collision_markers(capsules_world[closest_cap_i],
                                       capsules_world[closest_cap_j],
                                       best_result, *c_self->marker_pub);
    }

    return constraint_value;
}
// ============================================================================
// HELPER: CALCOLA JACOBIANO DI UN PUNTO SU UNA CAPSULA
// ============================================================================

Eigen::MatrixXd compute_capsule_jacobian(
    const Eigen::Matrix4d &T_link,
    const Capsule &cap_def,
    const Eigen::MatrixXd &J_link,
    double t_param,
    int NJ)
{
    // Calcola posizioni A e B in world frame
    Eigen::Matrix4d T = T_link * cap_def.T_offset;
    Eigen::Vector3d A_world = T.block<3, 1>(0, 3);
    Eigen::Vector3d B_world = A_world + T.block<3, 1>(0, 2) * cap_def.length;

    // Origine del link
    Eigen::Vector3d p_link = T_link.block<3, 1>(0, 3);

    // Jacobiani traslativi e rotativi del link
    Eigen::MatrixXd J_trans = J_link.block(0, 0, 3, NJ);
    Eigen::MatrixXd J_rot = J_link.block(3, 0, 3, NJ);

    // Bracci di leva
    Eigen::Vector3d r_A = A_world - p_link;
    Eigen::Vector3d r_B = B_world - p_link;

    // Jacobiani geometrici per A e B
    Eigen::MatrixXd J_A(3, NJ), J_B(3, NJ);
    for (int col = 0; col < NJ; col++)
    {
        Eigen::Vector3d omega = J_rot.col(col);
        J_A.col(col) = J_trans.col(col) + omega.cross(r_A);
        J_B.col(col) = J_trans.col(col) + omega.cross(r_B);
    }

    // Interpolazione lineare basata su t_param
    return (1.0 - t_param) * J_A + t_param * J_B;
}

// ============================================================================
// HELPER: VISUALIZZAZIONE SELF-COLLISION IN RVIZ
// ============================================================================

void publish_self_collision_markers(
    const CapsuleWorld &capA,
    const CapsuleWorld &capB,
    const CapsuleDistanceResult &result,
    ros::Publisher pub)
{
    visualization_msgs::MarkerArray marker_array;

    // Marker per capsula A (rossa)
    visualization_msgs::Marker markerA;
    markerA.header.frame_id = "panda_link0";
    markerA.header.stamp = ros::Time::now();
    markerA.ns = "self_collision_capsule_A";
    markerA.id = 0;
    markerA.type = visualization_msgs::Marker::CYLINDER;
    markerA.action = visualization_msgs::Marker::ADD;

    // Posizione: punto medio tra A e B della capsula
    markerA.pose.position.x = (capA.A.x() + capA.B.x()) / 2.0;
    markerA.pose.position.y = (capA.A.y() + capA.B.y()) / 2.0;
    markerA.pose.position.z = (capA.A.z() + capA.B.z()) / 2.0;
    markerA.pose.orientation.w = 1.0;

    markerA.scale.x = capA.radius * 2.0;
    markerA.scale.y = capA.radius * 2.0;
    markerA.scale.z = (capA.B - capA.A).norm();

    markerA.color.r = 1.0;
    markerA.color.g = 0.0;
    markerA.color.b = 0.0;
    markerA.color.a = 0.5;

    marker_array.markers.push_back(markerA);

    // Marker per capsula B (blu)
    visualization_msgs::Marker markerB = markerA;
    markerB.ns = "self_collision_capsule_B";
    markerB.id = 1;
    markerB.pose.position.x = (capB.A.x() + capB.B.x()) / 2.0;
    markerB.pose.position.y = (capB.A.y() + capB.B.y()) / 2.0;
    markerB.pose.position.z = (capB.A.z() + capB.B.z()) / 2.0;
    markerB.scale.x = capB.radius * 2.0;
    markerB.scale.y = capB.radius * 2.0;
    markerB.scale.z = (capB.B - capB.A).norm();
    markerB.color.r = 0.0;
    markerB.color.b = 1.0;

    marker_array.markers.push_back(markerB);

    // Linea tra punti più vicini
    visualization_msgs::Marker line;
    line.header = markerA.header;
    line.ns = "self_collision_distance";
    line.id = 2;
    line.type = visualization_msgs::Marker::LINE_STRIP;
    line.action = visualization_msgs::Marker::ADD;

    geometry_msgs::Point p1, p2;
    p1.x = result.p_capsule.x();
    p1.y = result.p_capsule.y();
    p1.z = result.p_capsule.z();

    p2.x = result.p_obstacle.x();
    p2.y = result.p_obstacle.y();
    p2.z = result.p_obstacle.z();

    line.points.push_back(p1);
    line.points.push_back(p2);

    line.scale.x = 0.005;
    line.color.r = 1.0;
    line.color.g = 1.0;
    line.color.b = 0.0;
    line.color.a = 1.0;

    marker_array.markers.push_back(line);

    pub.publish(marker_array);
}