#pragma once

#include <vector>
#include <Eigen/Dense>

// Assumiamo che la struct Capsule sia definita qui o inclusa da un altro file.
// Se l'hai già definita altrove, rimuovi questa definizione.
struct Capsule {
    int link_index;
    double radius;
    double length;
    Eigen::Matrix4d T_offset;
};

struct Fr3CollisionModel {

    static std::vector<Capsule> get_definitions() {
        std::vector<Capsule> capsule_definitions;

        // --- fr3_link0 (Index 0) ---
        {
            Capsule cap;
            cap.link_index = 0;
            cap.radius = 0.055000;
            cap.length = 0.030000;
            cap.T_offset << 0.0000, 0.0000, 1.0000, -0.0750,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            -1.0000, 0.0000, 0.0000, 0.0600,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link1 (Index 1) ---
        {
            Capsule cap;
            cap.link_index = 1;
            cap.radius = 0.060000;
            cap.length = 0.283000;
            cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, -0.1915,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link2 (Index 2) ---
        {
            Capsule cap;
            cap.link_index = 2;
            cap.radius = 0.070000;
            cap.length = 0.120000;
            cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, 0.0000,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link3 (Index 3) ---
        {
            Capsule cap;
            cap.link_index = 3;
            cap.radius = 0.090000;
            cap.length = 0.150000;
            cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, -0.1450,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link4 (Index 4) ---
        {
            Capsule cap;
            cap.link_index = 4;
            cap.radius = 0.090000;
            cap.length = 0.120000;
            cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, 0.0000,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link5 (Index 5) ---
        {
            Capsule cap;
            cap.link_index = 5;
            cap.radius = 0.090000;
            cap.length = 0.100000;
            cap.T_offset << 1.0000, 0.0000, 0.0000, 0.0000,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, -0.2600,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        {
            Capsule cap;
            cap.link_index = 5;
            cap.radius = 0.055000;
            cap.length = 0.140000;
            cap.T_offset << 0.9968, -0.0799, 0.0000, 0.0000,
                            0.0799, 0.9968, 0.0000, 0.0800,
                            0.0000, 0.0000, 1.0000, -0.1300,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link6 (Index 6) ---
        {
            Capsule cap;
            cap.link_index = 6;
            cap.radius = 0.080000;
            cap.length = 0.080000;
            // Offset X corretto a -0.0100 per centrare sul giunto
            cap.T_offset << 1.0000, 0.0000, 0.0000, -0.0100,
                            0.0000, 1.0000, 0.0000, 0.0000,
                            0.0000, 0.0000, 1.0000, -0.0100,
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        // --- fr3_link7 (Index 7)
        {
            // Parte 1: Braccio Orizzontale
            Capsule cap;
            cap.link_index = 7;
            cap.radius = 0.050000;  
            cap.length = 0.120000;  
                
            cap.T_offset << 0.0000, 0.0000, 1.0000, 0.0600,  
                            0.0000, 1.0000, 0.0000, 0.0000,  
                            -1.0000, 0.0000, 0.0000, 0.0100, 
                            0.0000, 0.0000, 0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        {
            // Parte 2: Flangia Verticale/Obliqua
            Capsule cap;
            cap.link_index = 7;
            cap.radius = 0.050000; 
            cap.length = 0.150000; 

            cap.T_offset << 1.0000, 0.0000,  0.0000, 0.1050, 
                            0.0000, 0.0000, -1.0000, 0.0000, 
                            0.0000, 1.0000,  0.0000, 0.0300,
                            0.0000, 0.0000,  0.0000, 1.0000;
            capsule_definitions.push_back(cap);
        }

        return capsule_definitions;
    }
};