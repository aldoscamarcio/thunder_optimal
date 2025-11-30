#pragma once

#include <memory>
#include <string>
#include <vector>

#include <controller_interface/multi_interface_controller.h>
#include <hardware_interface/joint_command_interface.h>
#include <hardware_interface/robot_hw.h>
#include <ros/node_handle.h>
#include <ros/time.h>
#include <sensor_msgs/JointState.h>
#include <franka_hw/franka_model_interface.h>
#include <franka_hw/franka_state_interface.h>

#include <Eigen/Dense>

// ================= ACADOS HEADERS =================
// Includiamo le definizioni di Acados. 
// Nota: Se il compilatore si lamenta di tipi sconosciuti (es. ocp_nlp_config), 
// potrebbe essere necessario includere anche "acados_c/ocp_nlp_interface.h" 
// o "acados_c/external_function_interface.h" se non sono inclusi automaticamente dal file generato.
extern "C" {
#include "acados_solver_frankino_pos_ctrl.h"
}
// ==================================================

namespace panda_controllers {

class AcadosTorqueController : public controller_interface::MultiInterfaceController<
                                           franka_hw::FrankaModelInterface,
                                           hardware_interface::EffortJointInterface,
                                           franka_hw::FrankaStateInterface> {
 public:
  // CORREZIONE: MultiInterfaceController richiede root_nh E controller_nh
  bool init(hardware_interface::RobotHW* robot_hw, ros::NodeHandle& root_nh, ros::NodeHandle& controller_nh) override;
  
  void starting(const ros::Time&) override;
  void update(const ros::Time&, const ros::Duration& period) override;
  void stopping(const ros::Time&) override;

 private:
  // Funzione ausiliaria per saturare le coppie
  Eigen::Matrix<double, 7, 1> saturateTorqueRate(
      const Eigen::Matrix<double, 7, 1>& tau_d_calculated,
      const Eigen::Matrix<double, 7, 1>& tau_J_d);

  // Callback per il target
  void setCommandCB(const sensor_msgs::JointStateConstPtr& msg);

  // Handle Franka
  std::unique_ptr<franka_hw::FrankaModelHandle> model_handle_;
  std::unique_ptr<franka_hw::FrankaStateHandle> state_handle_;
  std::vector<hardware_interface::JointHandle> joint_handles_;

  // ROS
  ros::NodeHandle cvc_nh;
  ros::Subscriber sub_command_;
  ros::Publisher pub_err_;

  // Variabili di stato robot
  Eigen::Matrix<double, 7, 1> q_curr;
  Eigen::Matrix<double, 7, 1> dot_q_curr;
  Eigen::Matrix<double, 7, 1> tau_J_d; // Coppia desiderata letta dal robot (senza gravità)

  // Variabili Target (Command)
  Eigen::Matrix<double, 7, 1> command_q_d;      // Target posizione

  // ACADOS DATA STRUCTURES
  frankino_pos_ctrl_solver_capsule* acados_ocp_capsule;
  ocp_nlp_config* nlp_config;
  ocp_nlp_dims* nlp_dims;
  ocp_nlp_in* nlp_in;
  ocp_nlp_out* nlp_out;
  ocp_nlp_solver* nlp_solver;
  void* nlp_opts;

  // Buffer per Acados
  double acados_x0[14];   // Stato iniziale [q, q_dot]
  double acados_u0[7];    // Input ottimo [tau]
  double acados_p[7];     // Parametri [q_target]
  
  // Limiti sicurezza e costanti
  const double kDeltaTauMax{1.0};
  Eigen::Matrix<double, 7, 1> tau_limit;
};

}  // namespace panda_controllers