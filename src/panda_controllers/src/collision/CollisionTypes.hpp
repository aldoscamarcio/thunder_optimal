#pragma once
#include <Eigen/Dense>

enum class ObstacleType { PLANE, RECTANGLE, BOX, CAPSULE };

// Capsula definita nel frame di un link
struct CapsuleLocal {
    int link_index;
    double radius;
    double length;
    Eigen::Matrix4d T_offset;
};

// Capsula espressa nel mondo
struct CapsuleWorld {
    Eigen::Vector3d A;
    Eigen::Vector3d B;
    double radius;
};

// Risultato distanza
struct CapsuleDistanceResult {
    double distance;
    Eigen::Vector3d p_capsule;
    Eigen::Vector3d p_obstacle;
    Eigen::Vector3d normal;
    double t_capsule;
    double t_obstacle;
};

// Definizione ostacoli
struct Plane {
    Eigen::Vector3d P0;
    Eigen::Vector3d n;
};

struct Rectangle {
    Eigen::Vector3d P0;
    Eigen::Vector3d Ux;
    Eigen::Vector3d Uy;
    double width;
    double height;
};

struct BoxOBB {
    Eigen::Vector3d center;
    Eigen::Matrix3d axes;
    Eigen::Vector3d half_ext;
};

struct Obstacle {
    ObstacleType type;
    Plane plane;
    Rectangle rect;
    BoxOBB box;
    CapsuleWorld capsule;
};
