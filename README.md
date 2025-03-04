# MATLAB Robotics Simulation

This repository contains MATLAB code for simulating flexible robotic arms with various configurations and solving their motion equations using numerical methods. The primary goal of the simulations is to analyze and visualize the deformation and movement of flexible rods under various forces.

## Files

### 1. `fig9.m`
This script is used to simulate and visualize the deformation of three flexible rods under specified forces and boundary conditions. It applies symbolic equations to model the system's behavior, solving for the rods' deflection and calculating various forces and moments acting on the system.[View the PDF](pdf/last1.pdf)



- **Parameters**: Tolerance for convergence, maximum iterations, Young's modulus, and geometry of the rods.
- **Key Outputs**: Deflection of rods, forces, and moments at different stages of deformation.

### 2. `fig11.m`
This script builds on `fig9.m` but incorporates additional complexities, such as symbolic force components and improved visualization techniques. It visualizes the deformation of the rods for various values of a parameter (`alpha`) and uses Newton's method to iteratively solve for the displacement.[View the PDF](pdf/last2.pdf)


- **Parameters**: Includes additional symbolic components for forces and moments.
- **Key Outputs**: Detailed plots of rod deformation and the forces at the ends of each rod.

### 3. `fig12.m`
This file further extends the previous models by adding more variables and improving the equation system to represent real-world robotic arm setups. The main difference is the inclusion of more rigid components and complex boundary conditions.
[View the PDF](pdf/last3.pdf)

- **Parameters**: Similar to the previous files, but adds more variables related to the force and moment at each joint.
- **Key Outputs**: Includes detailed results for multiple configurations and force distributions.

## How to Use

### Prerequisites
To run the scripts, you will need:
- MATLAB with **Symbolic Math Toolbox** enabled
- Basic knowledge of **MATLAB scripting**

### Running the Scripts
1. Download or clone this repository to your local machine.
2. Open MATLAB and navigate to the folder where the files are saved.
3. Run each script (`fig9.m`, `fig11.m`, or `fig12.m`) individually to simulate different robot arm configurations and visualize the results.

### Example Usage
```matlab
>> run('fig9.m');

