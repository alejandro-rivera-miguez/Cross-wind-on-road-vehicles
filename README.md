# Heavy Vehicle Dynamic Stability Under Crosswind 🌬️🚛

[![Vehicle Dynamics](https://img.shields.io/badge/Simulation-Dynamic_Model-blue.svg)]()
[![Control](https://img.shields.io/badge/Controller-PID_Driver_Model-orange.svg)]()
[![Academic](https://img.shields.io/badge/University-Politecnico%20di%20Milano-blue.svg)](https://www.polimi.it/)

A comprehensive parametric study evaluating the dynamic stability and rollover risk of heavy goods vehicles subjected to severe crosswind conditions. 

Developed for the **Wind Engineering** course at Politecnico di Milano (A.A. 2024-2025).

---

## 📌 Project Overview

Heavy goods vehicles are highly susceptible to strong lateral winds, particularly on exposed infrastructures like bridges and highways. Crosswind-induced rollovers and severe lane departures pose significant safety risks. 

This project develops a dynamic vehicle model integrating aerodynamic loading, payload variations, and a closed-loop driver model to predict critical rollover thresholds and assess lateral trajectory deviations.

### 🎯 Key Objectives & Scope
* **Aerodynamic Modeling:** Simulation of turbulent crosswind profiles with varying mean velocities (25 & 30 m/s) and turbulence intensities (7%, 14%, 25%).
* **Vehicle Dynamics:** Analysis of critical parameters including loading factors ($\lambda$) and roll stiffness distribution ($\tau$) between front and rear axles.
* **Driver Modeling:** Implementation of a PID controller to simulate driver steering corrections during wind gusts.
* **Maneuver Simulation:** Evaluating stability not just on straight paths, but during complex maneuvers such as single/double lane changes (overtaking) and cornering.

## 📂 Repository Contents

* **`docs/`**: Contains the final academic report and the slide deck detailing the mathematical models, aerodynamic assumptions, and conclusions.
* **`src/`**: Contains the simulation scripts and dynamic models used to perform the parametric sweeps and plot the stability boundaries.

## 📈 Key Findings
1. **Payload Dependency:** The rollover threshold exhibits a non-linear behavior; instability is highest at low payloads due to reduced tyre grip, while dynamic behavior improves as weight increases.
2. **Roll Stiffness:** Optimal rollover resistance is achieved at intermediate values of the front-to-total roll stiffness ratio, highlighting the importance of suspension tuning.
3. **Maneuvering Effects:** Steering inputs during turns and lane changes drastically reduce the critical rollover speed, proving that combined crosswind and path-following scenarios are the most critical for heavy vehicle safety.

## 👨‍💻 Authors (Group N)

* **Francisco Javier Martín López**
* **Alejandro Rivera Míguez**
* **Mikel Segovia Díaz**

---
Professor: **Prof. Daniele Rocchi, Prof. Federico Zanelli**  
Politecnico di Milano

---
*Disclaimer: This is an academic project. The results are derived from a simplified multi-body dynamic model intended for educational evaluation of wind engineering concepts.*
