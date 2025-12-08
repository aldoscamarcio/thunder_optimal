// DistanceFunctions.cpp
// Implementazioni ottimizzate per le funzioni di distanza capsula <-> primitiva
// Basato su: Christer Ericson - Real-Time Collision Detection (capitoli su segment-segment, distances)

#include "DistanceFunctions.hpp"
#include "CollisionTypes.hpp"

#include <cmath>
#include <array>
#include <limits>

// Costanti e utility per la gestione numerica
static constexpr double EPS_DBL = 1e-12;  // Tolleranza per confronti con zero
static constexpr double INF_DBL = std::numeric_limits<double>::infinity();  // Valore infinito

// Funzioni utility inline per ottimizzazione
static inline double sqr(double x) { return x * x; }  // Calcola il quadrato di un numero
static inline double clamp01(double v) { return (v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v); }  // Limita valore tra 0 e 1

// -----------------------------
// closestSegmentSegment (Ericson) - ottimizzata e robusta
// Calcola i punti più vicini tra due segmenti 3D e restituisce i parametri e le coordinate
// -----------------------------
void closestSegmentSegment(
    const Eigen::Vector3d &P0, const Eigen::Vector3d &P1,  // Endpoints primo segmento
    const Eigen::Vector3d &Q0, const Eigen::Vector3d &Q1,  // Endpoints secondo segmento
    double &s_out, double &t_out,                          // Parametri di output lungo i segmenti
    Eigen::Vector3d &ptP, Eigen::Vector3d &ptQ)            // Punti più vicini di output
{
    // Basato sull'algoritmo di Ericson per punti più vicini tra segmenti
    Eigen::Vector3d u = P1 - P0;  // Vettore direzione del segmento P
    Eigen::Vector3d v = Q1 - Q0;  // Vettore direzione del segmento Q
    Eigen::Vector3d w = Q0 - P0;  // Vettore che connette l'inizio di P all'inizio di Q (standard)
    
    // Calcolo dei prodotti scalari per il sistema lineare
    double a = u.dot(u);      // Quadrato della lunghezza di u (sempre >= 0)
    double b = u.dot(v);      // Proiezione di u su v (misura parallelismo tra segmenti)
    double c = v.dot(v);      // Quadrato della lunghezza di v (sempre >= 0)
    double d = u.dot(w);      // Proiezione di w su u (u·(Q0 - P0))
    double e = v.dot(w);      // Proiezione di w su v (v·(Q0 - P0))
    double D = a * c - b * b; // Determinante del sistema (denominatore)

    // Variabili per il calcolo dei parametri lungo i segmenti
    double sc, sN, sD = D; // Parametro s = sN / sD per segmento P
    double tc, tN, tD = D; // Parametro t = tN / tD per segmento Q

    // Calcolo dei parametri dei punti più vicini
    if (D < EPS_DBL)
    {             // Segmenti quasi paralleli
        sN = 0.0; // Forza l'uso di s = 0 sul segmento P
        sD = 1.0; // Evita divisione per zero
        tN = e;   // Parametro per t basato solo sul segmento Q
        tD = c;
    }
    else
    {
        // Calcolo dei numeratori usando la regola di Cramer
        // Formule corrette con w = Q0 - P0:
        sN = (c * d - b * e); // Numeratore per s: (c·d - b·e)
        tN = (a * e - b * d); // Numeratore per t: (a·e - b·d)

        // Clamp di sN nell'intervallo [0, sD] (assicura 0 ≤ s ≤ 1)
        if (sN < 0.0)
        {
            sN = 0.0; // Forza al punto iniziale di P (s = 0)
            tN = e;   // Ricalcola tN per s = 0
            tD = c;
        }
        else if (sN > sD)
        {
            sN = sD;    // Forza al punto finale di P (s = 1)
            tN = e + b; // Ricalcola tN per s = 1: v·(Q0 - P1) = v·(w - u) = e - b
            tD = c;
        }
    }

    // Clamp di tN nell'intervallo [0, tD] (assicura 0 ≤ t ≤ 1)
    if (tN < 0.0)
    {
        tN = 0.0; // Forza al punto iniziale di Q (t = 0)
        // Ricalcola sN per questo t (trovare punto su P più vicino a Q0)
        // Proiezione di (P0 - Q0) su u = -d (perché d = u·(Q0 - P0))
        if (-d < 0.0)
            sN = 0.0;           // Proiezione negativa → usa P0
        else if (-d > a)
            sN = sD;            // Proiezione > lunghezza → usa P1
        else
        {
            sN = -d;            // Proiezione tra 0 e a
            sD = a;             // Normalizza per ‖u‖²
        }
    }
    else if (tN > tD)
    {
        tN = tD; // Forza al punto finale di Q (t = 1)
        // Ricalcola sN per questo t (trovare punto su P più vicino a Q1)
        // Proiezione di (P0 - Q1) su u = u·(P0 - Q0 - v) = -d - b
        if ((-d - b) < 0.0)
            sN = 0.0;           // Proiezione negativa → usa P0
        else if ((-d - b) > a)
            sN = sD;            // Proiezione > lunghezza → usa P1
        else
        {
            sN = (-d - b);      // Proiezione tra 0 e a
            sD = a;             // Normalizza per ‖u‖²
        }
    }

    // Calcolo finale dei parametri con controllo di tolleranza numerica
    sc = (std::abs(sN) < EPS_DBL ? 0.0 : sN / sD);
    tc = (std::abs(tN) < EPS_DBL ? 0.0 : tN / tD);

    // Assegnazione degli output
    s_out = sc;  // Parametro s normalizzato (0 ≤ s ≤ 1)
    t_out = tc;  // Parametro t normalizzato (0 ≤ t ≤ 1)

    // Calcolo delle coordinate dei punti più vicini
    ptP = P0 + u * sc;  // Punto più vicino sul segmento P: P(s) = P0 + s*(P1-P0)
    ptQ = Q0 + v * tc;  // Punto più vicino sul segmento Q: Q(t) = Q0 + t*(Q1-Q0)
}

// -----------------------------
// dist_capsule_capsule
// Calcola la distanza signed tra due capsule usando l'algoritmo segmento-segmento
// Distanza signed: >0 separate, <0 penetranti, =0 si toccano
// -----------------------------
double dist_capsule_capsule(
    const CapsuleWorld &A,  // Prima capsula
    const CapsuleWorld &B,  // Seconda capsula
    CapsuleDistanceResult *out)  // Risultato opzionale dettagliato
{
    double s, t;  // Parametri lungo le capsule
    Eigen::Vector3d pA, pB;  // Punti più vicini
    // Trova i punti più vicini tra i segmenti centrali delle capsule
    closestSegmentSegment(A.A, A.B, B.A, B.B, s, t, pA, pB);

    // Calcolo della distanza tra i punti più vicini
    Eigen::Vector3d diff = pA - pB;
    double d = diff.norm();  // Distanza euclidea
    double radii = A.radius + B.radius;  // Somma dei raggi
    double signedDist = d - radii;  // Distanza signed

    if (out)
    {
        // Popolazione della struttura di output
        out->distance = signedDist;
        out->p_capsule = pA;  // Punto sulla prima capsula
        out->p_obstacle = pB; // Punto sulla seconda capsula
        out->t_capsule = s;   // Parametro lungo la prima capsula
        out->t_obstacle = t;  // Parametro lungo la seconda capsula
        
        // Calcolo della normale di collisione
        if (d > 1e-12)
            out->normal = diff / d;  // Normalizza il vettore differenza
        else
        {
            // Caso degenere: punti coincidenti, calcola normale alternativa
            Eigen::Vector3d u = A.B - A.A;  // Asse della capsula
            Eigen::Vector3d alt = u.cross(Eigen::Vector3d::UnitX());  // Prova con asse X
            if (alt.squaredNorm() < 1e-12)
                alt = u.cross(Eigen::Vector3d::UnitY());  // Prova con asse Y
            if (alt.squaredNorm() < 1e-12)
                out->normal = Eigen::Vector3d::UnitZ();  // Fallback a asse Z
            else
                out->normal = alt.normalized();  // Usa vettore perpendicolare normalizzato
        }
    }
    return signedDist;
}

// -----------------------------
// dist_capsule_plane 
// Calcola la distanza signed tra una capsula e un piano
// Gestisce sia il caso di penetrazione che di separazione
// -----------------------------
double dist_capsule_plane(
    const CapsuleWorld &cap,  // Capsula
    const Plane &pl,          // Piano
    CapsuleDistanceResult *out)  // Risultato opzionale
{
    // Calcolo distanze signed degli endpoints dal piano
    double dA = pl.n.dot(cap.A - pl.P0);
    double dB = pl.n.dot(cap.B - pl.P0);
    
    // CASO SPECIALE: Controlla se il segmento è parallelo al piano
    // (per evitare divisione per zero e gestire caso speciale)
    Eigen::Vector3d u = cap.B - cap.A;
    double u_dot_n = u.dot(pl.n);
    
    if (std::abs(u_dot_n) < EPS_DBL) {
        // Segmento parallelo al piano: tutti i punti hanno la stessa distanza
        double d = dA;  // dA = dB per segmenti paralleli
        double signedDist = std::abs(d) - cap.radius;  // Distanza signed
        
        if (out) {
            // Usa il punto centrale come rappresentativo
            Eigen::Vector3d center = (cap.A + cap.B) * 0.5;
            Eigen::Vector3d proj = center - pl.n * d;  // Proiezione sul piano
            
            out->distance = signedDist;
            out->p_capsule = center;     // Punto sulla capsula (centro)
            out->p_obstacle = proj;      // Proiezione sul piano
            out->t_capsule = 0.5;        // Parametro del centro
            out->t_obstacle = 0.0;
            // Normale punta dal piano verso la capsula
            out->normal = (d >= 0.0) ? pl.n : -pl.n;
        }
        return signedDist;
    }
    
    // CASO 1: La capsula INTERSECA il piano (dA e dB hanno segni opposti)
    if (dA * dB <= 0.0) {
        // La capsula attraversa il piano → punto più vicino è l'intersezione
        // La distanza minima tra asse e piano è 0, quindi considerando il raggio:
        double signedDist = -cap.radius;  // Penetrazione minima = raggio
        
        // Trova il punto di intersezione tra segmento e piano
        // t = distanza da A / lunghezza totale proiettata sulla normale
        double t = dA / (dA - dB);  // Parametro di intersezione
        t = clamp01(t);  // Assicura che sia tra 0 e 1 (per sicurezza numerica)
        Eigen::Vector3d P_intersect = cap.A + t * (cap.B - cap.A);  // Punto di intersezione
        
        if (out) {
            out->distance = signedDist;
            out->p_capsule = P_intersect;   // Punto sulla capsula (intersezione)
            out->p_obstacle = P_intersect;  // Punto sul piano (stessa posizione)
            out->t_capsule = t;             // Parametro di intersezione lungo capsula
            out->t_obstacle = 0.0;
            // Normale punta verso l'esterno rispetto alla capsula
            // (verso il lato dove si trova l'endpoint con distanza positiva)
            out->normal = (dA >= 0.0) ? pl.n : -pl.n;
        }
        return signedDist;
    }
    
    // CASO 2: A e B della capsula sono dallo stesso lato del piano
    else {
        // Scegli l'endpoint geometricamente più vicino al piano (in valore assoluto)
        bool useA = std::abs(dA) <= std::abs(dB);
        double d = useA ? dA : dB;           // Distanza signed del punto più vicino
        Eigen::Vector3d P = useA ? cap.A : cap.B;  // Punto più vicino sulla capsula
        Eigen::Vector3d proj = P - pl.n * d; // Proiezione ortogonale sul piano
        
        // Distanza signed = distanza minima tra asse e piano, meno il raggio
        double signedDist = std::abs(d) - cap.radius;
        
        if (out) {
            out->distance = signedDist;
            out->p_capsule = P;             // Punto sulla capsula (estremo più vicino)
            out->p_obstacle = proj;         // Sua proiezione ortogonale sul piano
            out->t_capsule = useA ? 0.0 : 1.0;  // 0 per A, 1 per B
            out->t_obstacle = 0.0;
            // Normale punta dal piano verso la capsula
            // (se d positivo: capsula sul lato positivo della normale)
            out->normal = (d >= 0.0) ? pl.n : -pl.n;
        }
        return signedDist;
    }
}

// -----------------------------
// dist_capsule_rectangle
// Calcola la distanza tra una capsula e un rettangolo 3D
// Considera sia la faccia che i bordi del rettangolo
// -----------------------------
double dist_capsule_rectangle(
    const CapsuleWorld &cap, 
    const Rectangle &R, 
    CapsuleDistanceResult *out) 
{
    const double EPSILON = 1e-12;
    double minD2 = INF_DBL; // Distanza minima al quadrato
    
    // Output temporanei per tracciare la feature più vicina
    Eigen::Vector3d closest_P_cap, closest_P_rect;
    double best_t_cap = 0.0;
    
    // 1. Calcolo normale e assi locali
    // Ericson Sez 5.1.4.2: Rectangle definito da centro e assi [cite: 261]
    // Qui manteniamo la tua struttura P0 + u*w + v*h
    Eigen::Vector3d n = R.Ux.cross(R.Uy);
    double n_norm = n.norm();
    if (n_norm < EPSILON) return -INF_DBL; // Degenerato
    n /= n_norm;

    // ---------------------------------------------------------
    // FASE A: Verifica Intersezione (Penetrazione profonda)
    // Ericson Sez 5.3.1: Intersect Segment Against Plane [cite: 1376]
    // Se la capsula attraversa la faccia, la distanza è 0 (o negativa)
    // ---------------------------------------------------------
    Eigen::Vector3d d_cap = cap.B - cap.A;
    double denom = n.dot(d_cap);
    
    if (std::abs(denom) > EPSILON) {
        double t = n.dot(R.P0 - cap.A) / denom;
        if (t >= 0.0 && t <= 1.0) {
            Eigen::Vector3d P_int = cap.A + t * d_cap;
            // Verifica se il punto di intersezione è dentro il rettangolo (Point in Polygon)
            Eigen::Vector3d v = P_int - R.P0;
            double u_proj = v.dot(R.Ux);
            double v_proj = v.dot(R.Uy);
            
            if (u_proj >= 0.0 && u_proj <= R.width && 
                v_proj >= 0.0 && v_proj <= R.height) {
                
                // Penetrazione rilevata
                if (out) {
                    out->distance = -cap.radius; // O profondità reale se desiderata
                    out->p_capsule = P_int;
                    out->p_obstacle = P_int;
                    out->normal = (denom < 0) ? n : -n; // Normale opposta alla direzione incidente
                    out->t_capsule = t;
                    out->t_obstacle = 0.0; // Generico per interno
                }
                return -cap.radius;
            }
        }
    }

    // ---------------------------------------------------------
    // FASE B: Estremi della Capsula contro Interno Faccia (Face Region)
    // Ericson Sez 5.1.10: Case (a) Endpoint vs Interior 
    // ---------------------------------------------------------
    const Eigen::Vector3d* cap_endpoints[2] = { &cap.A, &cap.B };
    double cap_t_vals[2] = { 0.0, 1.0 };

    for (int i = 0; i < 2; ++i) {
        const Eigen::Vector3d& P = *cap_endpoints[i];
        // Proietta P sul piano del rettangolo: Q = P - dist * n
        double dist_plane = n.dot(P - R.P0);
        Eigen::Vector3d Q = P - dist_plane * n;

        // Verifica se la proiezione Q è dentro i bordi del rettangolo
        Eigen::Vector3d v = Q - R.P0;
        double u_proj = v.dot(R.Ux);
        double v_proj = v.dot(R.Uy);

        if (u_proj >= 0.0 && u_proj <= R.width && 
            v_proj >= 0.0 && v_proj <= R.height) {
            
            double d2 = dist_plane * dist_plane;
            if (d2 < minD2) {
                minD2 = d2;
                closest_P_cap = P;
                closest_P_rect = Q;
                best_t_cap = cap_t_vals[i];
            }
        }
    }

    // ---------------------------------------------------------
    // FASE C: Segmento Capsula contro 4 Bordi (Edge Region)
    // Ericson Sez 5.1.9: Closest Points of Two Line Segments [cite: 619]
    // Ericson Sez 5.1.10: Case (b) Segment vs Edge [cite: 850]
    // ---------------------------------------------------------
    Eigen::Vector3d p00 = R.P0;
    Eigen::Vector3d p10 = R.P0 + R.Ux * R.width;
    Eigen::Vector3d p01 = R.P0 + R.Uy * R.height;
    Eigen::Vector3d p11 = p10 + R.Uy * R.height;

    const Eigen::Vector3d *edges[4][2] = {
        {&p00, &p10}, {&p10, &p11}, {&p11, &p01}, {&p01, &p00}
    };

    for (int i = 0; i < 4; ++i) {
        double s, t;
        Eigen::Vector3d pc, pr;
        // closestSegmentSegment gestisce parallelismo e clamping
        closestSegmentSegment(cap.A, cap.B, *edges[i][0], *edges[i][1], s, t, pc, pr);
        
        double d2 = (pc - pr).squaredNorm();
        if (d2 < minD2) {
            minD2 = d2;
            closest_P_cap = pc;
            closest_P_rect = pr;
            best_t_cap = s;
        }
    }

    // ---------------------------------------------------------
    // Calcolo Risultato Finale
    // ---------------------------------------------------------
    double dist = std::sqrt(minD2);
    double signedDist = dist - cap.radius;

    if (out) {
        out->distance = signedDist;
        out->p_capsule = closest_P_cap;
        out->p_obstacle = closest_P_rect;
        out->t_capsule = best_t_cap;
        // out->t_obstacle difficile da definire uniformemente (bordo vs faccia), lascio 0 o calcolo locale
        
        if (dist > EPSILON) {
            out->normal = (closest_P_cap - closest_P_rect) / dist;
        } else {
            // Caso di contatto o penetrazione leggera 
            // Fallback alla normale del piano
            out->normal = n; 
        }
    }

    return signedDist;
}

// -----------------------------
// dist_capsule_box (OBB) ottimizzata e robusta
// Strategia:
//  - Trasforma la capsula in coordinate locali della box
//  - Test veloce: endpoints clammati alla box (punto->AABB)
//  - Valuta segmento vs tutti i 12 bordi della box
//  - Restituisce la migliore distanza
// -----------------------------
double dist_capsule_box(
    const CapsuleWorld &cap,  // Capsula
    const BoxOBB &box,        // Box orientata (OBB)
    CapsuleDistanceResult *out)  // Risultato opzionale
{
    // Trasformazione: mondo -> coordinate locali della box
    Eigen::Matrix3d R = box.axes.transpose();  // Rotazione mondo->locale
    Eigen::Vector3d A_local = R * (cap.A - box.center);  // Trasforma endpoint A
    Eigen::Vector3d B_local = R * (cap.B - box.center);  // Trasforma endpoint B

    // Funzione helper per clammare punti all'interno della AABB locale
    auto clampLocal = [&](const Eigen::Vector3d &p)
    {
        return Eigen::Vector3d(
            std::clamp(p.x(), -box.half_ext.x(), box.half_ext.x()),  // Clamp X
            std::clamp(p.y(), -box.half_ext.y(), box.half_ext.y()),  // Clamp Y
            std::clamp(p.z(), -box.half_ext.z(), box.half_ext.z())); // Clamp Z
    };

    // Inizializzazione per la ricerca della distanza minima
    double best = INF_DBL;  // Miglior distanza trovata
    Eigen::Vector3d bestSegLocal = Eigen::Vector3d::Zero();  // Miglior punto sulla capsula (locale)
    Eigen::Vector3d bestBoxLocal = Eigen::Vector3d::Zero();  // Miglior punto sulla box (locale)
    double bestT = 0.0;  // Miglior parametro sulla capsula

    // 1) Test sugli endpoints della capsula
    {
        Eigen::Vector3d pts[2] = {A_local, B_local};  // Endpoints in coordinate locali
        for (int i = 0; i < 2; ++i)
        {
            const Eigen::Vector3d &P = pts[i];
            Eigen::Vector3d Q = clampLocal(P);  // Punto più vicino sulla box
            double dcur = (P - Q).norm();       // Distanza endpoint-box
            if (dcur < best)
            {
                best = dcur;
                bestSegLocal = P;
                bestBoxLocal = Q;
                bestT = (i == 0) ? 0.0 : 1.0;  // 0 per A, 1 per B
            }
        }
    }

    // 2) Test su tutti i bordi della box (12 bordi totali)
    Eigen::Vector3d he = box.half_ext;  // Half-extents per comodità
    
    // Vertici locali della box (tutte le combinazioni di ±half-extents)
    std::array<Eigen::Vector3d, 8> v;
    v[0] = Eigen::Vector3d(-he.x(), -he.y(), -he.z());  // Vertice posteriore-inferiore-sinistro
    v[1] = Eigen::Vector3d( he.x(), -he.y(), -he.z());  // Vertice posteriore-inferiore-destro
    v[2] = Eigen::Vector3d( he.x(),  he.y(), -he.z());  // Vertice posteriore-superiore-destro
    v[3] = Eigen::Vector3d(-he.x(),  he.y(), -he.z());  // Vertice posteriore-superiore-sinistro
    v[4] = Eigen::Vector3d(-he.x(), -he.y(),  he.z());  // Vertice anteriore-inferiore-sinistro
    v[5] = Eigen::Vector3d( he.x(), -he.y(),  he.z());  // Vertice anteriore-inferiore-destro
    v[6] = Eigen::Vector3d( he.x(),  he.y(),  he.z());  // Vertice anteriore-superiore-destro
    v[7] = Eigen::Vector3d(-he.x(),  he.y(),  he.z());  // Vertice anteriore-superiore-sinistro

    // Definizione degli edge della box come coppie di vertici
    const int E[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Faccia posteriore (4 bordi)
        {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Faccia anteriore (4 bordi)
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Bordi verticali (4 bordi)
    };

    // Per ogni bordo della box, calcola la distanza minima dalla capsula
    for (int ei = 0; ei < 12; ++ei)
    {
        Eigen::Vector3d e0 = v[E[ei][0]];  // Primo vertice del bordo
        Eigen::Vector3d e1 = v[E[ei][1]];  // Secondo vertice del bordo
        double s_local, t_local;
        Eigen::Vector3d segPtLocal, edgePtLocal;
        // Trova punti più vicini tra capsula e bordo corrente
        closestSegmentSegment(A_local, B_local, e0, e1, s_local, t_local, segPtLocal, edgePtLocal);
        double dcur = (segPtLocal - edgePtLocal).norm();  // Distanza calcolata
        if (dcur < best)
        {
            best = dcur;
            bestSegLocal = segPtLocal;
            bestBoxLocal = edgePtLocal;
            bestT = s_local;  // Parametro sulla capsula
        }
    }

    // Trasformazione inversa: coordinate locali -> mondo
    Eigen::Vector3d segPtWorld = box.axes * bestSegLocal + box.center;  // Punto capsula in mondo
    Eigen::Vector3d boxPtWorld = box.axes * bestBoxLocal + box.center;  // Punto box in mondo

    double signedDist = best - cap.radius;  // Distanza signed considerando il raggio

    if (out) {
        out->distance = signedDist;
        out->t_capsule = bestT;
        out->p_capsule = segPtWorld;
        out->p_obstacle = boxPtWorld;
        
        // Calcolo della normale di collisione
        Eigen::Vector3d diff = segPtWorld - boxPtWorld;
        double dn = diff.norm();
        if (dn > 1e-12) {
            out->normal = diff / dn;  // Normalizza la differenza
        } else {
            out->normal = Eigen::Vector3d::UnitX();  // Fallback a asse X
        }

        out->t_obstacle = 0.0;  // Non utilizzato per le box
    }

    return signedDist;
}

// -----------------------------
// Dispatcher principale
// Instrada alla funzione di distanza appropriata in base al tipo di ostacolo
// -----------------------------
double capsule_distance(const CapsuleWorld &cap, const Obstacle &obs, CapsuleDistanceResult *out)
{
    switch (obs.type)
    {
    case ObstacleType::PLANE:
        return dist_capsule_plane(cap, obs.plane, out);      // Capsula vs Piano
    case ObstacleType::RECTANGLE:
        return dist_capsule_rectangle(cap, obs.rect, out);   // Capsula vs Rettangolo
    case ObstacleType::BOX:
        return dist_capsule_box(cap, obs.box, out);          // Capsula vs Box
    case ObstacleType::CAPSULE:
        return dist_capsule_capsule(cap, obs.capsule, out);  // Capsula vs Capsula
    default:
        return INF_DBL;  // Tipo di ostacolo sconosciuto
    }
}