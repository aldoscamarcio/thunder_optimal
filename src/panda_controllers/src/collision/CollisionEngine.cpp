#include "CollisionEngine.hpp"

double CollisionEngine::distance(
    const CapsuleWorld &cap,
    const Obstacle &obs,
    CapsuleDistanceResult *res)
{
    return capsule_distance(cap, obs, res);
}

Eigen::RowVectorXd CollisionEngine::gradient(
    const CapsuleLocal &cap_local,
    const Eigen::Matrix4d &T_link,
    const Eigen::MatrixXd &J_link,
    const CapsuleDistanceResult &res,
    int NJ)
{
    return capsule_distance_gradient(cap_local, T_link, J_link, res, NJ);
}
