#include <iostream>
#include <Eigen/Dense>

#include "collision/CollisionTypes.hpp"
#include "collision/CapsuleModels.hpp"
#include "collision/DistanceFunctions.hpp"

void printResult(const std::string &name, const CapsuleDistanceResult &r)
{
    std::cout << "==== " << name << " ====\n";
    std::cout << "Distance:       " << r.distance << "\n";
    std::cout << "Capsule point:  " << r.p_capsule.transpose() << "\n";
    std::cout << "Obstacle point: " << r.p_obstacle.transpose() << "\n";
    std::cout << "Normal:         " << r.normal.transpose() << "\n";
    std::cout << "t_capsule:      " << r.t_capsule << "\n";
    std::cout << "=========================\n\n";
}
int main()
{
    std::cout << "===== COLLISION VALIDATION TEST =====\n\n";

    //--------------------------------------------------------------------------------------------------
    // TEST 1 — Capsule → Plane
    // Piano z = 0, normale +Z
    // Capsula con estremi A=(0,0,1), B=(0,0,2), radius=0.2
    // Distanza attesa: |1| - 0.2 = 0.8
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap1;
    cap1.A = Eigen::Vector3d(0, 0, 1);
    cap1.B = Eigen::Vector3d(0, 0, 2);
    cap1.radius = 0.2;

    Plane plane;
    plane.P0 = Eigen::Vector3d(0, 0, 0);
    plane.n = Eigen::Vector3d(0, 0, 1);

    CapsuleDistanceResult r1;
    dist_capsule_plane(cap1, plane, &r1);
    printResult("Capsule -> Plane", r1);

    //--------------------------------------------------------------------------------------------------
    // TEST 2 — Capsule → Rectangle (Inside)
    // Rettangolo nel piano z=0, di dimensioni 2x2 centrato in (1,1,0)
    // Capsula con A=(1,1,1), B=(1,1,3), radius=0.3
    // Distanza attesa: |1| - 0.3 = 0.7
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap2;
    cap2.A = Eigen::Vector3d(1, 1, 1);
    cap2.B = Eigen::Vector3d(1, 1, 3);
    cap2.radius = 0.3;

    Rectangle rect2;
    rect2.P0 = Eigen::Vector3d(0, 0, 0); // angolo basso sinistro
    rect2.Ux = Eigen::Vector3d(1, 0, 0);
    rect2.Uy = Eigen::Vector3d(0, 1, 0);
    rect2.width = 2.0;
    rect2.height = 2.0;

    CapsuleDistanceResult r2;
    dist_capsule_rectangle(cap2, rect2, &r2);
    printResult("Capsule -> Rectangle (Inside)", r2);

    //--------------------------------------------------------------------------------------------------
    // TEST 3 — Capsule → OBB (Box)
    // Box centrato in (2,0,0), dimensioni 2x2x2
    // Capsula A=(0,0,0), B=(0,0,1), radius=0.2
    // Distanza attesa: distanza capsula->box approx = (2 - 0.2) = 1.8
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap3;
    cap3.A = Eigen::Vector3d(0, 0, 0);
    cap3.B = Eigen::Vector3d(0, 0, 1);
    cap3.radius = 0.2;

    BoxOBB box;
    box.center = Eigen::Vector3d(2, 0, 0);
    box.axes = Eigen::Matrix3d::Identity();  // orientazione identità
    box.half_ext = Eigen::Vector3d(1, 1, 1); // box 2x2x2

    CapsuleDistanceResult r3;
    dist_capsule_box(cap3, box, &r3);
    printResult("Capsule -> Box (OBB)", r3);

    //--------------------------------------------------------------------------------------------------
    // TEST 4 — Capsule → Capsule
    // Capsule 1: A=(0,0,0), B=(1,0,0), r=0.2
    // Capsule 2: A=(2,0,0), B=(3,0,0), r=0.2
    // Distanza attesa: center-to-center = 1, distance = 1 - (0.2+0.2) = 0.6
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap4_A;
    cap4_A.A = Eigen::Vector3d(0, 0, 0);
    cap4_A.B = Eigen::Vector3d(1, 0, 0);
    cap4_A.radius = 0.2;

    CapsuleWorld cap4_B;
    cap4_B.A = Eigen::Vector3d(2, 0, 0);
    cap4_B.B = Eigen::Vector3d(3, 0, 0);
    cap4_B.radius = 0.2;

    CapsuleDistanceResult r4;
    dist_capsule_capsule(cap4_A, cap4_B, &r4);
    printResult("Capsule -> Capsule", r4);

    //--------------------------------------------------------------------------------------------------
    // TEST 5 — Capsule → Rectangle (Outside Corner)
    // Rettangolo P0=(0,0,0), 2x2
    // Capsula A=(3,3,1), B=(3,3,2), r=0.2
    // Distanza attesa: Proiezione (3,3,0) -> Vicino (2,2,0)
    // Dist = sqrt((3-2)^2 + (3-2)^2 + (1-0)^2) - 0.2 = sqrt(3) - 0.2 = 1.732 - 0.2 = 1.532 (CORRETTO)
    // Nota: Nel commento originale c'era un calcolo diverso, qui usiamo la logica standard Euclidea
    // Punto capsula più vicino al rettangolo è A(3,3,1).
    // Punto rettangolo più vicino ad A è (2,2,0).
    // Distanza centri = sqrt(1^2 + 1^2 + 1^2) = sqrt(3) = ~1.732
    // Distanza superficie = 1.732 - 0.2 = 1.532
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap5;
    cap5.A = Eigen::Vector3d(3, 3, 1); // sopra l'angolo superiore destro, esterno
    cap5.B = Eigen::Vector3d(3, 3, 2);
    cap5.radius = 0.2;

    Rectangle rect5;
    rect5.P0 = Eigen::Vector3d(0, 0, 0);
    rect5.Ux = Eigen::Vector3d(1, 0, 0);
    rect5.Uy = Eigen::Vector3d(0, 1, 0);
    rect5.width = 2.0;
    rect5.height = 2.0;

    CapsuleDistanceResult r5;
    dist_capsule_rectangle(cap5, rect5, &r5);
    printResult("Capsule -> Rectangle (Outside Corner)", r5);

    //--------------------------------------------------------------------------------------------------
    // TEST 6 — Capsule → Rectangle (Outside Side)
    // Rettangolo P0=(0,0,0), 2x2
    // Capsula A=(3,1,1), B=(3,1,2), r=0.1
    // Distanza attesa:
    // Proiezione su piano: (3,1). Punto rettangolo più vicino (2,1). Dist orizzontale = 1.
    // Altezza z=1.
    // Distanza centri = sqrt(1^2 + 1^2) = sqrt(2) = 1.414
    // Distanza superficie = 1.414 - 0.1 = 1.314
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap6;
    cap6.A = Eigen::Vector3d(3, 1, 1); // fuori dal lato destro
    cap6.B = Eigen::Vector3d(3, 1, 2);
    cap6.radius = 0.1;

    Rectangle rect6;
    rect6.P0 = Eigen::Vector3d(0, 0, 0);
    rect6.Ux = Eigen::Vector3d(1, 0, 0);
    rect6.Uy = Eigen::Vector3d(0, 1, 0);
    rect6.width = 2.0;
    rect6.height = 2.0;

    CapsuleDistanceResult r6;
    dist_capsule_rectangle(cap6, rect6, &r6);
    printResult("Capsule -> Rectangle (Outside Side)", r6);

    //--------------------------------------------------------------------------------------------------
    // TEST 7 — Capsule → Rectangle (Above & Outside)
    // Rettangolo P0=(0,0,0), 2x2
    // Capsula A=(3,3,5), B=(3,3,6), r=0.5
    // Distanza attesa:
    // Proiezione (3,3). Punto rettangolo vicino (2,2). Distanza planare sqrt(2).
    // Altezza z=5.
    // Distanza centri = sqrt(1^2 + 1^2 + 5^2) = sqrt(27) = 5.196
    // Distanza superficie = 5.196 - 0.5 = 4.696
    //--------------------------------------------------------------------------------------------------
    CapsuleWorld cap7;
    cap7.A = Eigen::Vector3d(3, 3, 5);
    cap7.B = Eigen::Vector3d(3, 3, 6);
    cap7.radius = 0.5;

    Rectangle rect7;
    rect7.P0 = Eigen::Vector3d(0, 0, 0);
    rect7.Ux = Eigen::Vector3d(1, 0, 0);
    rect7.Uy = Eigen::Vector3d(0, 1, 0);
    rect7.width = 2.0;
    rect7.height = 2.0;

    CapsuleDistanceResult r7;
    dist_capsule_rectangle(cap7, rect7, &r7);
    printResult("Capsule -> Rectangle (Above & Outside)", r7);


    //--------------------------------------------------------------------------------------------------
    // TEST 8 — Capsule → Rectangle (Oblique & Outside)
    // Rettangolo P0=(0,0,0), 2x2 (quindi lato destro a x=2, y=[0,2])
    // Capsula "incrociata": fissa a x=3, ma inclinata in Z.
    // A=(3, 0, 2), B=(3, 2, 0). Raggio = 0.2.
    //
    // Analisi:
    // La capsula passa "a fianco" del rettangolo.
    // Il punto più vicino del rettangolo è il centro del lato destro: (2, 1, 0).
    // Il punto più vicino della capsula (per simmetria) è il suo punto medio: (3, 1, 1).
    // Distanza Euclidea centri = sqrt((3-2)^2 + (1-1)^2 + (1-0)^2) 
    //                          = sqrt(1 + 0 + 1) = sqrt(2) ≈ 1.4142
    // Distanza Superficie attesa = 1.4142 - 0.2 = 1.2142
    //--------------------------------------------------------------------------------------------------
    
    CapsuleWorld cap8;
    cap8.A = Eigen::Vector3d(3, 0, 2); 
    cap8.B = Eigen::Vector3d(3, 2, 0);
    cap8.radius = 0.2;

    Rectangle rect8;
    rect8.P0 = Eigen::Vector3d(0, 0, 0);
    rect8.Ux = Eigen::Vector3d(1, 0, 0);
    rect8.Uy = Eigen::Vector3d(0, 1, 0);
    rect8.width = 2.0;
    rect8.height = 2.0;

    CapsuleDistanceResult r8;
    dist_capsule_rectangle(cap8, rect8, &r8);
    printResult("Capsule -> Rectangle (Oblique & Outside)", r8);

    //--------------------------------------------------------------------------------------------------
    // TEST 9 — Capsule → Rectangle (Diagonal Cross over Corner)
    // Rettangolo 2x2 a z=0. Angolo in alto a destra è (2,2,0).
    // Capsula sospesa a z=1, che "taglia" l'angolo a 45 gradi.
    // A=(3, 1, 1), B=(1, 3, 1). Raggio = 0.2.
    //
    // Analisi:
    // Punto Medio Capsula M = (2, 2, 1).
    // Distanza M -> Angolo Rettangolo (2,2,0) è esattamente 1.0 (lungo Z).
    //
    // Verifica Estremi:
    // A=(3,1,1). Punto rettangolo più vicino è (2,1,0). Dist = sqrt(1^2 + 0 + 1^2) = 1.414.
    // B=(1,3,1). Punto rettangolo più vicino è (1,2,0). Dist = sqrt(0 + 1^2 + 1^2) = 1.414.
    //
    // Poiché 1.0 < 1.414, il punto più vicino DEVE essere il centro della capsula (t=0.5).
    // Distanza attesa = 1.0 - 0.2 = 0.8.
    //--------------------------------------------------------------------------------------------------

    CapsuleWorld cap9;
    cap9.A = Eigen::Vector3d(3, 1, 1);
    cap9.B = Eigen::Vector3d(1, 3, 1);
    cap9.radius = 0.2;

    Rectangle rect9;
    rect9.P0 = Eigen::Vector3d(0, 0, 0);
    rect9.Ux = Eigen::Vector3d(1, 0, 0);
    rect9.Uy = Eigen::Vector3d(0, 1, 0);
    rect9.width = 2.0;
    rect9.height = 2.0;

    CapsuleDistanceResult r9;
    dist_capsule_rectangle(cap9, rect9, &r9);
    printResult("Capsule -> Rectangle (Diagonal Cross over Corner)", r9);

    //--------------------------------------------------------------------------------------------------
    // TEST 10 — Capsule → Box (Edge Interaction)
    // Box 2x2x2 centrato in origine. Spigolo superiore destro a (x=1, z=1).
    // Capsula orizzontale parallela allo spigolo, spostata lungo Y.
    // A=(2, 0, 1), B=(2, 2, 1). Raggio 0.2.
    //
    // Analisi:
    // Il Box va da x=[-1,1], z=[-1,1].
    // La capsula è a x=2, z=1. Altezza costante z=1 (livello cima box).
    // Distanza orizzontale tra la linea della capsula (x=2) e lo spigolo del box (x=1) è 1.0.
    // Il punto più vicino sul box è lungo lo spigolo (1, y, 1).
    // Distanza Euclidea = 1.0.
    // Distanza Superficie = 1.0 - 0.2 = 0.8.
    //--------------------------------------------------------------------------------------------------

    CapsuleWorld cap10;
    cap10.A = Eigen::Vector3d(2, 0, 1);
    cap10.B = Eigen::Vector3d(2, 2, 1);
    cap10.radius = 0.2;

    BoxOBB box10;
    box10.center = Eigen::Vector3d(0, 0, 0);
    box10.half_ext = Eigen::Vector3d(1, 1, 1);
    box10.axes = Eigen::Matrix3d::Identity();

    CapsuleDistanceResult r10;
    dist_capsule_box(cap10, box10, &r10);
    printResult("Capsule -> Box (Edge Interaction)", r10);

    //--------------------------------------------------------------------------------------------------
    // TEST 11 — Capsule → Rotated Box (The "Diamond" setup)
    // Box 2x2x2 ruotato di 45 gradi su Z.
    // Visto dall'alto è un rombo. Il vertice più a destra è a X = sqrt(2) ≈ 1.414.
    // Capsula verticale ferma a X=2.0, Y=0.
    //
    // Analisi:
    // Semiextent = 1. La distanza dal centro all'angolo (nel piano XY) è sqrt(1^2 + 1^2) = 1.4142.
    // La capsula è a x=2.0.
    // Distanza dal centro capsula all'angolo del box = 2.0 - 1.4142 = 0.5858.
    // Distanza superficie = 0.5858 - 0.2 (raggio) = 0.3858.
    //--------------------------------------------------------------------------------------------------

    CapsuleWorld cap11;
    cap11.A = Eigen::Vector3d(2, 0, -1);
    cap11.B = Eigen::Vector3d(2, 0, 1);
    cap11.radius = 0.2;

    BoxOBB box11;
    box11.center = Eigen::Vector3d(0, 0, 0);
    box11.half_ext = Eigen::Vector3d(1, 1, 1);
    // Rotazione 45 gradi su Z
    box11.axes = Eigen::AngleAxisd(M_PI / 4.0, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    CapsuleDistanceResult r11;
    dist_capsule_box(cap11, box11, &r11);
    printResult("Capsule -> Rotated Box (Diamond)", r11);

    //--------------------------------------------------------------------------------------------------
    // TEST 12 — Capsule → Box Face (Oblique Projection)
    // Box standard. Faccia superiore a z=1.
    // Capsula inclinata sospesa sopra la faccia superiore.
    // A=(0, 0, 2), B=(0, 2, 3). Raggio 0.2.
    //
    // Analisi:
    // Il punto A(0,0,2) è direttamente sopra il centro della faccia superiore (0,0,1).
    // La distanza verticale è 1.0.
    // Il punto B(0,2,3) è sopra lo spigolo (0,1,1) ma molto più in alto (z=3 vs z=1 -> dist 2).
    // Il punto più vicino è inequivocabilmente A.
    // Distanza = 1.0 - 0.2 = 0.8.
    // Questo verifica che l'algoritmo proietti correttamente sulla faccia e non cerchi spigoli errati.
    //--------------------------------------------------------------------------------------------------

    CapsuleWorld cap12;
    cap12.A = Eigen::Vector3d(0, 0, 2);
    cap12.B = Eigen::Vector3d(0, 2, 3);
    cap12.radius = 0.2;

    BoxOBB box12;
    box12.center = Eigen::Vector3d(0, 0, 0);
    box12.half_ext = Eigen::Vector3d(1, 1, 1);
    box12.axes = Eigen::Matrix3d::Identity();

    CapsuleDistanceResult r12;
    dist_capsule_box(cap12, box12, &r12);
    printResult("Capsule -> Box Face (Oblique Projection)", r12);

    std::cout << "===== TEST COMPLETATI =====\n";
    return 0;
}