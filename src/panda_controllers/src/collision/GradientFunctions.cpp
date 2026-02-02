#include "GradientFunctions.hpp"

Eigen::RowVectorXd capsule_distance_gradient(
    const CapsuleLocal &cap_local,
    const Eigen::Matrix4d &T_link,
    const Eigen::MatrixXd &J_link,
    const CapsuleDistanceResult &res,
    int NJ)
{
    /*******************************************************
     * 1) Ricostruisci i punti della capsula in world frame
     *******************************************************/
    // 1) Ottieni la trasformazione globale della capsula
        Eigen::Matrix4d T = T_link * cap_local.T_offset;

        // 2) M è il centro della capsula (posizione), quindi un Vector3d (3x1)
        // T.block<3,1>(0,3) estrae le prime 3 righe della 4a colonna
        Eigen::Vector3d M = T.block<3, 1>(0, 2);

        // 3) Estrai l'asse lungo cui si sviluppa la capsula (solitamente l'asse Z, colonna 2)
        Eigen::Vector3d direction = T.block<3, 1>(0, 2);

        // 4) Calcola i punti estremi A e B
        // Sottraiamo/aggiungiamo metà lunghezza lungo la direzione dell'asse
        Eigen::Vector3d A = M - direction * (cap_local.length / 2.0);
        Eigen::Vector3d B = M + direction * (cap_local.length / 2.0);

        // Punto nel frame del link (se ti serve l'asse Z del link)
        Eigen::Vector3d p_link_z = T_link.block<3, 1>(0, 2);

    /*******************************************************
     * 2) Separa Jacobiano traslazionale e rotazionale
     *    J_link è 6×N = [ v ; ω ]
     *******************************************************/
    Eigen::MatrixXd Jv = J_link.topRows(3);    // 3 × NJ
    Eigen::MatrixXd Jw = J_link.bottomRows(3); // 3 × NJ

    /*******************************************************
     * 3) Costruisci Jacobiani del punto A e del punto B
     *
     *    p = p_link + R * offset
     *    J = Jv + ω × (p - p_link)
     *******************************************************/
    Eigen::Vector3d rA = A - p_link_z;
    Eigen::Vector3d rB = B - p_link_z;

    Eigen::MatrixXd J_A(3, NJ), J_B(3, NJ);

    for (int j = 0; j < NJ; j++)
    {
        Eigen::Vector3d w = Jw.col(j);

        J_A.col(j) = Jv.col(j) + w.cross(rA);
        J_B.col(j) = Jv.col(j) + w.cross(rB);
    }

    /*******************************************************
     * 4) Jacobiano del punto della capsula più vicino
     *    t_capsule ∈ [0,1] dal risultato distanza
     *******************************************************/
    Eigen::MatrixXd J_closest =
        (1.0 - res.t_capsule) * J_A +
         res.t_capsule        * J_B;

    /*******************************************************
     * 5) Gradiente della distanza rispetto a q:
     *
     *    ∂d/∂q = normalᵀ · J_closest
     *******************************************************/
    Eigen::RowVectorXd ddist_dq = res.normal.transpose() * J_closest;

    return ddist_dq;
}
