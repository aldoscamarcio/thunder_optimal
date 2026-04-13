# Nonlinear Trajectory Optimization for Robotic Manipulators

This repository contains the code and algorithms developed for my Master's Thesis in Robotics and Automation Engineering at the University of Pisa.

The project proposes an innovative framework based on Nonlinear Optimization for generating optimal trajectories for robotic manipulators in the presence of obstacles. Unlike traditional probabilistic planners (such as RRT* or PRM*) that often generate sub-optimal, purely geometric paths, this approach intrinsically integrates the physics and geometry of the system. It minimizes energy consumption while ensuring strict compliance with complex kinematic and dynamic constraints.


## 🎯 Key Features

    Arbitrary Boundary Conditions: Allows setting non-zero final velocities and accelerations, enabling the execution of complex dynamic tasks (e.g., throwing objects).

    Closed-Form Collision Modeling: Evaluates the minimum distance between the robot and the environment by approximating mechanical links with capsules and obstacles with convex primitives (spheres, cylinders, rectangles, planes).

    High Repeatability: Provides highly repeatable solutions under identical initial conditions, overcoming the strong stochastic dispersion typical of sampling-based methods.

## 🛠 Project Tasks

This repository reflects the following implementation phases described in the thesis:

- **Optimal Control Formulation:** Formulated the trajectory planning as a finite-horizon problem, defining joint accelerations as the optimization variables.- **Initialization (Minimum Jerk):** Developed a module to generate a minimum-jerk polynomial trajectory to be used as the initial guess for the optimizer.
- **Constraint Integration:**
  - *Kinematic and Dynamic:* Imposed physical limits on position, velocity, acceleration, and motor torques.
  - *Obstacle Avoidance:* Implemented an obstacle repulsion system by calculating analytical distances between volumes.
- **Simulation Validation:** Designed benchmark testing scenarios with increasing complexity: spherical obstacle, rectangular obstacle, and double obstacle.
- **Experimental Validation (Hardware):**
  - Bypassing a rectangular obstacle placed on the direct trajectory.
  - Executing a dynamic task of throwing a ball into a container, utilizing the qb SoftHand robotic hand and imposing a specific non-zero final velocity.
- **Statistical Comparative Analysis:** Acquired data at 1 kHz to measure Control Energy, Jerk, and Path Length. The data was statistically analyzed against RRT* and PRM* using Raincloud Plots and the Wilcoxon non-parametric test.

<video src="assets/MVI_8147_light.MP4" controls="controls" style="max-width: 730px;">
</video>

## 💻 Tools Used

- **NLopt (LD-MMA Solver):** Used to minimize the cost functional based on gradients (Method of Moving Asymptotes).
- **ROS (Robot Operating System):** Framework used to encapsulate instantaneous references into messages and for real-time manipulator control.
- **RViz:** Used for rendering and validating the approximated model (capsules) of the robot and obstacles in the workspace.
- **Hardware:** Robotic manipulator (Franka Emika Panda) equipped with a qb SoftHand end-effector.

## ⚙️ System Requirements

- **OS:** Ubuntu Linux 20.04.
- **Framework:** ROS (Robot Operating System) properly installed and configured.
- **Math Libraries:** NLopt library for solving nonlinear optimization problems.
- **Hardware Dependencies (Optional):** Specific ROS packages for manipulator control and qb SoftHand drivers (only required for execution on the physical robot).

  
  ### Author: Aldo Scamarcio
