#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdlib.h>

// Struttura dati
struct TrajectoryPoint {
    std::vector<double> q;
    std::vector<double> dq;
    std::vector<double> ddq;
};

bool loadTrajectory(const std::string& filename, std::vector<TrajectoryPoint>& trajectory) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        ROS_ERROR_STREAM("Errore: Impossibile aprire il file " << filename);
        return false;
    }

    std::string line;
    
    // --- MODIFICA IMPORTANTE: SALTARE L'HEADER ---
    // Leggiamo la prima riga (q1,q2,...) e la ignoriamo
    if (std::getline(file, line)) {
        ROS_INFO("Header rilevato e saltato: %s...", line.substr(0, 20).c_str());
    }
    // ---------------------------------------------

    int row_count = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value_str;
        std::vector<double> row_values;

        // Parsing separato da virgola (',')
        while (std::getline(ss, value_str, ',')) { 
            try {
                // Rimuove eventuali spazi bianchi residui
                if (!value_str.empty()) {
                    row_values.push_back(std::stod(value_str));
                }
            } catch (...) {
                ROS_WARN("Valore non numerico trovato alla riga %d, colonna ignorata.", row_count + 2);
            }
        }

        // Controllo: Il tuo CSV ha 21 colonne (7q + 7dq + 7ddq)
        if (row_values.size() >= 21) {
            TrajectoryPoint pt;
            
            // Assegnazione vettori
            pt.q.assign(row_values.begin(), row_values.begin() + 7);
            pt.dq.assign(row_values.begin() + 7, row_values.begin() + 14);
            pt.ddq.assign(row_values.begin() + 14, row_values.begin() + 21);
            
            trajectory.push_back(pt);
            row_count++;
        } else {
            if (!row_values.empty()) // Evita di spammare errori per righe vuote a fine file
                ROS_WARN("Riga %d scartata: trovati %lu valori invece di 21.", row_count + 2, row_values.size());
        }
    }
    
    file.close();
    ROS_INFO_STREAM("Caricati con successo " << row_count << " punti traiettoria.");
    return row_count > 0;
}

int main(int argc, char **argv) {
    ros::init(argc, argv, "csv_trajectory_player");
    ros::NodeHandle nh("~");

    // Publisher verso il topic del tuo Computed Torque Controller
    ros::Publisher pub_cmd = nh.advertise<sensor_msgs::JointState>("/computed_torque_controller/command", 1000);

    // Parametri dal launch file
    std::string csv_path;
    double frequency;
    
    nh.param<std::string>("csv_path", csv_path, ""); 
    nh.param<double>("frequency", frequency, 1000.0); 

    if (csv_path.empty()) {
        ROS_ERROR("Parametro csv_path non settato!");
        return -1;
    }

    ROS_INFO("Apro il file: %s", csv_path.c_str());

    std::vector<TrajectoryPoint> trajectory;
    if (!loadTrajectory(csv_path, trajectory)) {
        ROS_ERROR("Traiettoria vuota o file non valido. Esco.");
        return -1;
    }

    ROS_INFO("TRAIETTORIA PRONTA. Premi INVIO nel terminale per eseguire...");
    std::cin.ignore(); 

    ros::Rate loop_rate(frequency);
    sensor_msgs::JointState msg;
    
    size_t i = 0;
    while (ros::ok() && i < trajectory.size()) {
        
        msg.header.stamp = ros::Time::now();
        
        // Assegnazione q, dq, ddq (effort)
        msg.position = trajectory[i].q;
        msg.velocity = trajectory[i].dq;
        msg.effort   = trajectory[i].ddq; 

        pub_cmd.publish(msg);

        ros::spinOnce();
        loop_rate.sleep();
        i++;
    }

    ROS_INFO("Esecuzione terminata.");
    
    // Stop finale (opzionale): invia comando a velocità 0
    if (!trajectory.empty()) {
        msg.header.stamp = ros::Time::now();
        msg.position = trajectory.back().q;     // Resta nell'ultima posizione
        msg.velocity = std::vector<double>(7, 0.0); // Ferma velocità
        msg.effort   = std::vector<double>(7, 0.0); // Ferma accelerazione
        pub_cmd.publish(msg);
    }

    return 0;
}