#include "panda_controllers/acados_controller.h" 
#include <pluginlib/class_list_macros.h>
#include <algorithm>

namespace panda_controllers {

bool AcadosTorqueController::init(hardware_interface::RobotHW* robot_hw, ros::NodeHandle& node_handle) {
    this->cvc_nh = node_handle;

    // ===========================================================================
    // 1. PARAMETRI ROS (arm_id, joint_names)
    // ===========================================================================
    std::string arm_id;
    if (!node_handle.getParam("arm_id", arm_id)) {
        ROS_ERROR("AcadosTorqueController: Could not get parameter arm_id!");
        return false;
    }

    std::vector<std::string> joint_names;
    if (!node_handle.getParam("joint_names", joint_names) || joint_names.size() != 7) {
        ROS_ERROR("AcadosTorqueController: Error or Invalid size for joint_names!");
        return false;
    }

    // ===========================================================================
    // 2. INTERFACCIA MODELLO (FrankaModelInterface)
    // ===========================================================================
    franka_hw::FrankaModelInterface* model_interface = robot_hw->get<franka_hw::FrankaModelInterface>();
    if (model_interface == nullptr) {
        ROS_ERROR_STREAM("AcadosTorqueController: Error getting model interface from hardware!");
        return false;
    }

    try {
        model_handle_.reset(new franka_hw::FrankaModelHandle(model_interface->getHandle(arm_id + "_model")));
    } catch (hardware_interface::HardwareInterfaceException& ex) {
        ROS_ERROR_STREAM("AcadosTorqueController: Exception getting model handle from interface: " << ex.what());
        return false;
    }

    // ===========================================================================
    // 3. INTERFACCIA DI STATO (FrankaStateInterface)
    // ===========================================================================
    franka_hw::FrankaStateInterface* state_interface = robot_hw->get<franka_hw::FrankaStateInterface>();
    if (state_interface == nullptr) {
        ROS_ERROR_STREAM("AcadosTorqueController: Error getting state interface from hardware");
        return false;
    }

    try {
        state_handle_.reset(new franka_hw::FrankaStateHandle(state_interface->getHandle(arm_id + "_robot")));
    } catch (hardware_interface::HardwareInterfaceException& ex) {
        ROS_ERROR_STREAM("AcadosTorqueController: Exception getting state handle from interface: " << ex.what());
        return false;
    }

    // ===========================================================================
    // 4. INTERFACCIA DI SFORZO (EffortJointInterface)
    // ===========================================================================
    hardware_interface::EffortJointInterface* effort_joint_interface = robot_hw->get<hardware_interface::EffortJointInterface>();
    if (effort_joint_interface == nullptr) {
        ROS_ERROR_STREAM("AcadosTorqueController: Error getting effort joint interface from hardware!");
        return false;
    }

    // Puliamo il vettore per sicurezza
    joint_handles_.clear(); 
    
    for (size_t i = 0; i < 7; ++i) {
        try {
            joint_handles_.push_back(effort_joint_interface->getHandle(joint_names[i]));
        } catch (const hardware_interface::HardwareInterfaceException& ex) {
            ROS_ERROR_STREAM("AcadosTorqueController: Exception getting joint handles: " << ex.what());
            return false;
        }
    }

    // ===========================================================================
    // 5. INIZIALIZZAZIONE ACADOS
    // ===========================================================================
    ROS_INFO("Inizializzazione Solver Acados...");
    
    // Alloca la memoria per la capsula
    acados_ocp_capsule = frankino_pos_ctrl_acados_create_capsule();
    if (acados_ocp_capsule == nullptr) {
        ROS_ERROR("Impossibile creare la capsula Acados!");
        return false;
    }

    // Crea il solver
    int status = frankino_pos_ctrl_acados_create(acados_ocp_capsule);
    if (status) {
        ROS_ERROR_STREAM("Creazione solver Acados fallita con status: " << status);
        return false;
    }

    // Recupera i puntatori alle strutture interne
    nlp_config = frankino_pos_ctrl_acados_get_nlp_config(acados_ocp_capsule);
    nlp_dims = frankino_pos_ctrl_acados_get_nlp_dims(acados_ocp_capsule);
    nlp_in = frankino_pos_ctrl_acados_get_nlp_in(acados_ocp_capsule);
    nlp_out = frankino_pos_ctrl_acados_get_nlp_out(acados_ocp_capsule);
    nlp_solver = frankino_pos_ctrl_acados_get_nlp_solver(acados_ocp_capsule);
    nlp_opts = frankino_pos_ctrl_acados_get_nlp_opts(acados_ocp_capsule);

    // Inizializza ROS subscribers
    this->sub_command_ = node_handle.subscribe<sensor_msgs::JointState>("command", 1, &AcadosTorqueController::setCommandCB, this);
    this->pub_err_ = node_handle.advertise<sensor_msgs::JointState>("tracking_error", 1);

    // Default target
    command_q_d.setZero();
    
    ROS_INFO("AcadosController Inizializzato con successo.");
    return true;
}

void AcadosTorqueController::starting(const ros::Time& time) {
    // Ottieni stato iniziale del robot
    franka::RobotState robot_state = state_handle_->getRobotState();
    q_curr = Eigen::Map<Eigen::Matrix<double, 7, 1>>(robot_state.q.data());
    dot_q_curr = Eigen::Map<Eigen::Matrix<double, 7, 1>>(robot_state.dq.data());

    // Imposta il target iniziale sulla posizione corrente per evitare scatti
    command_q_d = q_curr;

    // Resetta il solver
    frankino_pos_ctrl_acados_reset(acados_ocp_capsule, 1);

    // Inizializza lo stato x0 e la guess iniziale
    for (int i = 0; i < 7; i++) {
        acados_x0[i] = q_curr[i];
        acados_x0[i + 7] = dot_q_curr[i];
    }

    // Imposta la guess iniziale (x) per tutti i nodi all'attuale stato stazionario
    for (int i = 0; i <= nlp_dims->N; i++) {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", acados_x0);
        
        // Imposta anche il parametro iniziale p (target)
        for(int j=0; j<7; j++) acados_p[j] = command_q_d[j];
        frankino_pos_ctrl_acados_update_params(acados_ocp_capsule, i, acados_p, 7);
    }
}

void AcadosTorqueController::update(const ros::Time&, const ros::Duration& period) {
    // 1. Leggi Stato Robot
    franka::RobotState robot_state = state_handle_->getRobotState();
    q_curr = Eigen::Map<Eigen::Matrix<double, 7, 1>>(robot_state.q.data());
    dot_q_curr = Eigen::Map<Eigen::Matrix<double, 7, 1>>(robot_state.dq.data());
    tau_J_d = Eigen::Map<Eigen::Matrix<double, 7, 1>>(robot_state.tau_J_d.data());

    // Prepara array x0 [q, q_dot]
    for (int i = 0; i < 7; i++) {
        acados_x0[i] = q_curr[i];
        acados_x0[i + 7] = dot_q_curr[i];
    }

    // 2. Imposta Vincoli Stato Iniziale (x0)
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", acados_x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", acados_x0);

    // 3. Imposta Parametri (Target)
    for (int i = 0; i < 7; i++) acados_p[i] = command_q_d[i];

    for (int i = 0; i <= nlp_dims->N; i++) {
        frankino_pos_ctrl_acados_update_params(acados_ocp_capsule, i, acados_p, 7);
    }

    // 4. Risolvi OCP
    int status = frankino_pos_ctrl_acados_solve(acados_ocp_capsule);

    if (status != ACADOS_SUCCESS) {
        ROS_WARN_THROTTLE(1.0, "Acados solver failed with status %d", status);
        // In caso di fallimento, potresti voler inviare coppie zero (gravità compensata) o mantenere l'ultimo comando valido.
        // Qui assumiamo che il robot gestisca la sicurezza o vada in gravity compensation.
        // return; 
    } else {
        // 5. Estrai Input Ottimo (u0)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", acados_u0);
    }

    Eigen::Matrix<double, 7, 1> tau_cmd_vector;
    for(int i=0; i<7; i++) tau_cmd_vector[i] = acados_u0[i];

    // 6. Saturazione e Sicurezza
    tau_cmd_vector = saturateTorqueRate(tau_cmd_vector, tau_J_d);

    // 7. Invia ai motori
    for (size_t i = 0; i < 7; i++) {
        joint_handles_[i].setCommand(tau_cmd_vector[i]);
    }

    // Debug / Logging Errori
    Eigen::Matrix<double, 7, 1> error = command_q_d - q_curr;
    Eigen::Matrix<double, 7, 1> dot_error = -dot_q_curr; 
    
    sensor_msgs::JointState error_msg;
    std::vector<double> err_vec(error.data(), error.data() + error.size());
    std::vector<double> dot_err_vec(dot_error.data(), dot_error.data() + dot_error.size());
    error_msg.header.stamp = ros::Time::now();
    error_msg.position = err_vec;
    error_msg.velocity = dot_err_vec;
    this->pub_err_.publish(error_msg);
}

void AcadosTorqueController::stopping(const ros::Time&) {
    // Gestione cleanup se necessaria
}

Eigen::Matrix<double, 7, 1> AcadosTorqueController::saturateTorqueRate(
    const Eigen::Matrix<double, 7, 1>& tau_d_calculated,
    const Eigen::Matrix<double, 7, 1>& tau_J_d) 
{
    Eigen::Matrix<double, 7, 1> tau_d_saturated = {};
    for (size_t i = 0; i < 7; i++) {
        double difference = tau_d_calculated[i] - tau_J_d[i];
        tau_d_saturated[i] = tau_J_d[i] + std::max(std::min(difference, kDeltaTauMax), -kDeltaTauMax);
    }
    return tau_d_saturated;
}

void AcadosTorqueController::setCommandCB(const sensor_msgs::JointStateConstPtr& msg) {
    if (msg->position.size() != 7) {
        ROS_ERROR_THROTTLE(1.0, "Target position size mismatch!");
        return;
    }
    command_q_d = Eigen::Map<const Eigen::Matrix<double, 7, 1>>(msg->position.data());
}

} // namespace panda_controllers

PLUGINLIB_EXPORT_CLASS(panda_controllers::AcadosTorqueController, controller_interface::ControllerBase);