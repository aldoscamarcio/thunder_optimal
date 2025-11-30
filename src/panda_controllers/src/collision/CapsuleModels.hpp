#pragma once
#include "CollisionTypes.hpp"

CapsuleWorld computeCapsuleWorld(
    const CapsuleLocal &cap,
    const Eigen::Matrix4d &T_link);
