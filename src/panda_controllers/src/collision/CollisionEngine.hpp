#pragma once
#include "CollisionTypes.hpp"
#include "CapsuleModels.hpp"
#include "DistanceFunctions.hpp"
#include "GradientFunctions.hpp"

class CollisionEngine {
public:
    CollisionEngine() = default;

    // Distanza capsula–ostacolo
    double distance(const CapsuleWorld &cap,
                    const Obstacle &obs,
                    CapsuleDistanceResult *res=nullptr);

    // Calcolo gradiente vincolo
    Eigen::RowVectorXd gradient(const CapsuleLocal &cap_local,
                                const Eigen::Matrix4d &T_link,
                                const Eigen::MatrixXd &J_link,
                                const CapsuleDistanceResult &res,
                                int NJ);
};
