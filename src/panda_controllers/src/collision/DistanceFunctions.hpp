#pragma once
#include "CollisionTypes.hpp"

// Utilità segment-segment
void closestSegmentSegment(
    const Eigen::Vector3d &P0, const Eigen::Vector3d &P1,
    const Eigen::Vector3d &Q0, const Eigen::Vector3d &Q1,
    double &s_out, double &t_out,
    Eigen::Vector3d &ptP, Eigen::Vector3d &ptQ);

// Distanze principali
double dist_capsule_plane(const CapsuleWorld &cap,
                          const Plane &pl,
                          CapsuleDistanceResult *out);

double dist_capsule_rectangle(const CapsuleWorld &cap,
                              const Rectangle &R,
                              CapsuleDistanceResult *out);

double dist_capsule_box(const CapsuleWorld &cap,
                        const BoxOBB &box,
                        CapsuleDistanceResult *out);

double dist_capsule_capsule(const CapsuleWorld &A,
                            const CapsuleWorld &B,
                            CapsuleDistanceResult *out);

// Dispatcher unico
double capsule_distance(const CapsuleWorld &cap,
                        const Obstacle &obs,
                        CapsuleDistanceResult *out = nullptr);
