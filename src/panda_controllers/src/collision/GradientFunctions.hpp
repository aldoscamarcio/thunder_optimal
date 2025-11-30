#pragma once
#include "CollisionTypes.hpp"
#include <Eigen/Dense>

Eigen::RowVectorXd capsule_distance_gradient(
    const CapsuleLocal &cap_local,
    const Eigen::Matrix4d &T_link,
    const Eigen::MatrixXd &J_link,  // 6×NJ
    const CapsuleDistanceResult &res,
    int NJ);
