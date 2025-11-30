#include "CapsuleModels.hpp"

CapsuleWorld computeCapsuleWorld(
    const CapsuleLocal &cap,
    const Eigen::Matrix4d &T_link)
{
    CapsuleWorld out;

    Eigen::Matrix4d T = T_link * cap.T_offset;

    out.A = T.block<3,1>(0,3);
    out.B = out.A + T.block<3,1>(0,2) * cap.length;
    out.radius = cap.radius;

    return out;
}
