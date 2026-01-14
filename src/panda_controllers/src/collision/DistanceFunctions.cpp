// DistanceFunctions.cpp
// Implementazioni ottimizzate per le funzioni di distanza capsula <-> primitiva
// Basato su: Christer Ericson - Real-Time Collision Detection (capitoli su segment-segment, distances)

#include "DistanceFunctions.hpp"
#include "CollisionTypes.hpp"
#include <cmath>
#include <array>
#include <limits>
#include <vector>

// Costanti e utility per la gestione numerica
static constexpr double EPS_DBL = 1e-12;  // Tolleranza per confronti con zero
static constexpr double INF_DBL = std::numeric_limits<double>::infinity();  // Valore infinito

// Funzioni utility inline per ottimizzazione
static inline double sqr(double x) { return x * x; }  // Calcola il quadrato di un numero
static inline double clamp01(double v) { return (v < 0.0) ? 0.0 : ((v > 1.0) ? 1.0 : v); }  // Limita valore tra 0 e 1

// Template generico per clamp
template<typename T>
static inline T clamp(T v, T minVal, T maxVal) { return (v < minVal) ? minVal : ((v > maxVal) ? maxVal : v); }

// -----------------------------
// closestSegmentSegment (Ericson) - ottimizzata e robusta
// Calcola i punti più vicini tra due segmenti 3D e restituisce i parametri e le coordinate
// -----------------------------

void closestSegmentSegment(
    const Eigen::Vector3d &P0, const Eigen::Vector3d &P1,
    const Eigen::Vector3d &Q0, const Eigen::Vector3d &Q1,
    double &s_out, double &t_out,
    Eigen::Vector3d &ptP, Eigen::Vector3d &ptQ)
{
    Eigen::Vector3d u = P1 - P0;
    Eigen::Vector3d v = Q1 - Q0;
    Eigen::Vector3d w = P0 - Q0;  // ✅ ATTENZIONE: P0 - Q0, non Q0 - P0
    
    double a = u.dot(u);
    double b = u.dot(v);
    double c = v.dot(v);
    double d = u.dot(w);
    double e = v.dot(w);
    double D = a * c - b * b;

    double sc, sN, sD = D;
    double tc, tN, tD = D;

    // Calcola i parametri non clampati
    if (D < EPS_DBL) {
        // Paralleli
        sN = 0.0;
        sD = 1.0;
        tN = e;
        tD = c;
    }
    else {
        // Caso generale
        sN = (b * e - c * d);
        tN = (a * e - b * d);

        // Clamp s
        if (sN < 0.0) {
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {
            sN = sD;
            tN = e + b;  // ✅ Con w = P0-Q0, questa è la formula corretta
            tD = c;
        }
    }

    // Clamp t
    if (tN < 0.0) {
        tN = 0.0;
        // Ricalcola s per t=0
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {
        tN = tD;
        // Ricalcola s per t=1
        if ((-d - b) < 0.0)
            sN = 0.0;
        else if ((-d - b) > a)
            sN = sD;
        else {
            sN = (-d - b);
            sD = a;
        }
    }

    sc = (std::abs(sN) < EPS_DBL ? 0.0 : sN / sD);
    tc = (std::abs(tN) < EPS_DBL ? 0.0 : tN / tD);

    s_out = sc;
    t_out = tc;

    ptP = P0 + u * sc;
    ptQ = Q0 + v * tc;
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
    const CapsuleWorld &cap,    // Capsula
    const Rectangle &R,         // Rettangolo
    CapsuleDistanceResult *out) // Risultato opzionale
{
    const double EPSILON = 1e-12;

    // 1. Calcolo della normale e setup assi
    Eigen::Vector3d n = R.Ux.cross(R.Uy).normalized();
    Eigen::Vector3d u = cap.B - cap.A; // Vettore asse capsula
    double u_dot_n = u.dot(n);

    // -----------------------------------------------------------
    // FASE A: Verifica collisione con la FACCIA (Face Region)
    // -----------------------------------------------------------
    
    // Trova il punto sulla linea della capsula più vicino al piano infinito del rettangolo
    double t_cap = 0.0;
    if (std::abs(u_dot_n) > EPSILON) {
        // Intersezione linea-piano
        t_cap = -n.dot(cap.A - R.P0) / u_dot_n;
        // Clamp tra 0 e 1 (per restare sul segmento della capsula)
        t_cap = std::max(0.0, std::min(1.0, t_cap));
    } else {
        // Parallelo: prendiamo il centro per stabilità
        t_cap = 0.5;
    }

    Eigen::Vector3d P_cap = cap.A + t_cap * u;     // Punto sulla capsula
    double d_plane = n.dot(P_cap - R.P0);          // Distanza signed dal piano
    Eigen::Vector3d P_proj = P_cap - n * d_plane;  // Punto proiettato sul piano

    // Calcolo coordinate locali 2D sul rettangolo
    Eigen::Vector3d v = P_proj - R.P0;
    double x_proj = v.dot(R.Ux);
    double y_proj = v.dot(R.Uy);

    // SE IL PUNTO CADE DENTRO I CONFINI DEL RETTANGOLO:
    if (x_proj >= 0.0 && x_proj <= R.width && y_proj >= 0.0 && y_proj <= R.height)
    {
        double distPlaneAbs = std::abs(d_plane);
        double signedDist = distPlaneAbs - cap.radius; // Distanza euclidea reale

        // Se siamo vicini alla collisione o dentro, attiviamo la logica di scivolamento
        if (out)
        {
            out->p_capsule = P_cap;
            out->p_obstacle = P_proj;
            out->t_capsule = t_cap;
            out->t_obstacle = 0.0; // Convenzionale per la faccia

            // --- LOGICA DI SCIVOLAMENTO (ANTI-BLOCCO) ---
            if (signedDist < 0) // Solo se c'è compenetrazione
            {
                // 1. Calcola distanza dai 4 bordi
                double d_left = x_proj;
                double d_right = R.width - x_proj;
                double d_bottom = y_proj;
                double d_top = R.height - y_proj;
                
                // Trova il bordo più vicino
                double min_edge_dist = std::min({d_left, d_right, d_bottom, d_top});
                
                // Direzione vettoriale verso quel bordo
                Eigen::Vector3d dir_edge = Eigen::Vector3d::Zero();
                if (min_edge_dist == d_left)       dir_edge = -R.Ux;
                else if (min_edge_dist == d_right) dir_edge = R.Ux;
                else if (min_edge_dist == d_bottom) dir_edge = -R.Uy;
                else                               dir_edge = R.Uy;

                // Fattore di pendenza laterale (0.5 - 0.8 è un buon range)
                double k_lateral = 0.7; 

                // FIX 1: Modifica il VALORE della distanza
                // Sottraiamo una penalità basata sulla distanza dal bordo.
                // Al centro del rettangolo il valore sarà molto più negativo che ai bordi.
                out->distance = signedDist - (min_edge_dist * k_lateral);

                // FIX 2: Modifica la NORMALE (Gradiente)
                // Sommiamo la normale del piano con la direzione verso l'uscita
                Eigen::Vector3d plane_n = (d_plane >= 0.0) ? n : -n;
                out->normal = (plane_n + k_lateral * dir_edge).normalized();
            }
            else 
            {
                // Nessuna collisione: Comportamento standard Euclideo
                out->distance = signedDist;
                out->normal = (d_plane >= 0.0) ? n : -n;
            }
        }
        
        // Ritorniamo il valore modificato se siamo in collisione (così l'ottimizzatore lo vede)
        return (out && signedDist < 0) ? out->distance : signedDist;
    }

    // -----------------------------------------------------------
    // FASE B: Collisione con i BORDI (Edge Region)
    // Se siamo qui, P_proj è fuori dal rettangolo 2D. 
    // Dobbiamo cercare la distanza minima segmento-segmento sui 4 lati.
    // -----------------------------------------------------------
    
    double minD2 = 1e18; // Infinito
    Eigen::Vector3d best_pc, best_pr;
    double best_s = 0.0, best_t = 0.0;
    
    // Vertici del rettangolo
    Eigen::Vector3d p00 = R.P0;
    Eigen::Vector3d p10 = R.P0 + R.Ux * R.width;
    Eigen::Vector3d p01 = R.P0 + R.Uy * R.height;
    Eigen::Vector3d p11 = p10 + R.Uy * R.height;

    // Array dei 4 bordi
    const Eigen::Vector3d *edges[4][2] = {
        {&p00, &p10}, // Bottom
        {&p10, &p11}, // Right
        {&p11, &p01}, // Top
        {&p01, &p00}  // Left
    };

    for (int i = 0; i < 4; ++i) {
        double s, t;
        Eigen::Vector3d pc, pr;
        // Funzione standard per distanza segmento-segmento (Assumendo tu l'abbia)
        closestSegmentSegment(cap.A, cap.B, *edges[i][0], *edges[i][1], s, t, pc, pr);

        double d2 = (pc - pr).squaredNorm();
        if (d2 < minD2) {
            minD2 = d2;
            best_pc = pc;
            best_pr = pr;
            best_s = s;
            best_t = t; // Parametro locale del bordo
        }
    }

    double distReal = std::sqrt(minD2);
    double signedDist = distReal - cap.radius;

    if (out) {
        out->distance = signedDist;
        out->p_capsule = best_pc;
        out->p_obstacle = best_pr;
        out->t_capsule = best_s;
        out->t_obstacle = best_t;
        
        // Normale standard per i bordi
        if (distReal > EPSILON) {
            out->normal = (best_pc - best_pr) / distReal;
        } else {
            // Caso degenere (tocco esatto): usiamo normale piano o una media
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
    Eigen::Vector3d A_local = box.axes.transpose() * (cap.A - box.center);
    Eigen::Vector3d B_local = box.axes.transpose() * (cap.B - box.center);
    
    // Funzione helper per clammare punti all'interno della AABB locale
    auto clampLocal = [&](const Eigen::Vector3d &p)
    {
        return Eigen::Vector3d(
            clamp(p.x(), -box.half_ext.x(), box.half_ext.x()),  // Clamp X
            clamp(p.y(), -box.half_ext.y(), box.half_ext.y()),  // Clamp Y
            clamp(p.z(), -box.half_ext.z(), box.half_ext.z())); // Clamp Z
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