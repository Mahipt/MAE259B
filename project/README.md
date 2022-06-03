# MAE259B __Team 8__ Project #
## Softrobot Motion Simulation in 2D Vascular-mimicking Network ##
- Abstract 
	- Soft continuum robots with active steerable probes capable of navigating in a complex vascular network hold great promise for medical applications. This repo contains a simulation of a soft robot controlled by an external magnetic field to navigate in a 2D vascular network. We established a DER model for the soft robot and successfully controlled the equivalent force exerted by the magnetic field to follow the desired trajectory. The next step of this project includes integration of DER mode and PD controller, soft robot model refinement, and determining collision between the robot and walls.  
	- Demo: 
		- ![picture alt](https://github.com/normanzmz/MAE259_Team8_Project/blob/main/Midterm_Progress/vascular_demo.jpg) 

## Table of Contents ## 
- General Info. 
- Techonologies
- Installation 
- Result
- Collaboration
- Reference

## 
### General Info: ###
- Class: __MAE259B__ (Advanced Topics in Solid Mechanics) 
- Team: __Team8__
- Author: 
	- Mingzhang Zhu 
	- Zhengqi Zhong 
	- YaoHsing Tseng

##
### Techonologies 
- DER Methods: 
- PD Control: 
- Force Computation: 
	- Viscosious Force: 
	- Contact Force with Maze: 
- Maze: 
	- Maze Generator: 
	- Maze contact function: 

##
### Installation 
- Environment requirements 
	- Matlab (version: 2022a) 
		- Simulink
	- Operating System: 
		- Windows 
		- Ubuntu 
		- MacOS
	
- Files and directories
	- Final/src
		- /main.m: Run this file to initialize all the parameters, and start the computation.
		- /Main_PD_Control.m: File to start PD controller simulation.  
		- /Control_Block.slx: (Simulink) 
			- Robot Systme (DER)
			- Control Loop 
			- Viscosious Force Computation. 
		- /maze.m: Create radom maze given the size of it. 
		- /load_maze.m: Load the existed maze file. 
		- /largemaze_fixed.mat: Generated maze (fixed, not random) 
		- /BoundaryJudge.m: Function to check whether nodes contact with the maze boundary (inside maze.)  
		- lib/*: Contains all computational libraries (written by Prof. M.Khalid Jawed.) 
- Install and Run: 
	1. Add Final/src/lib/* to Matlab path. 
	2. Run Final/src/main.m
	3. Run Final/src/Control_Block.slx 
	4. Simulation would start.  

##
###
- Result
	- Softrobot reacts to force in maze environment: 
		- ![picture alt](https://github.com/normanzmz/MAE259_Team8_Project/blob/main/Final/result/Contact%20Force%202.png)
	- Simulation video: 
		- https://user-images.githubusercontent.com/64833251/171540306-3c4d0a07-86de-4a3f-856c-9c63b8bb754c.mp4




##
### Collaboration 
- Prof. M.Khalid Jawed (Supervisor): 
	- Libraries
	- Guidance
- Authors: 
	- Mingzhang Zhu (Team Leader):   
		- Softrobot research and design.
		- Presentation
		- Force Computation: 
			- Viscosious Force
		- DER simulator design
	- Zhengqi Zhong: 
		- PD Controller Design. 
		- Robot Maze Contact Force Design. 
		- Simulation Environment Design. 
		- Presentation
		- Flow chart
		- Force Computation: 
			- Cotact Force. 
	- YaoHsing Tseng: 
		- Flow chart. 
		- Discussion
		- Github Readme Editing
		- Presentation.  


##
### Reference
1. Zion Tsz Ho Tse et al., Soft Robotics in Medical Applications, Journal of Medical Robotics Research, doi: 10.1142/S2424905X18410064
2. B. Vucelic et al., "The Aer-O-Scope Proof of Concept of a Pneumatic Skill-Independent, Self-Propelling Self- Navigating Colonoscope," Gastroenterology, vol. 130, no. 3, pp. 672-677, 2006.
3. H. G. Ren, X. & Tan, K. L., "Human-Compliant Body-Attached Soft Robots Towards Automatic Cooperative Ultrasound Imaging," 2016 20th IEEE International Conference on Computer Supported Cooperative Work in Design (CSCWD 2016), 2016.
4. Mark Runciman, Ara Darzi, and George P. Mylonas.Soft Robotics.Aug 2019.423-443. http://doi.org/10.1089/soro.2018.0136
5. M. Zhu, Y. Shen, A. J. Chiluisa, J. Song, L. Fichera and Y. Liu, "Optical Fiber Coupling System for Steerable Endoscopic Instruments," 2021 43rd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC), 2021, pp. 4871-4874, doi: 10.1109/EMBC46164.2021.9629658.
6. I. A. Chan, J. F. d’Almeida, A. J. Chiluisa, T. L. Carroll, Y. Liu, and L. Fichera, “On the merits of using angled fiber tips in office-based laser surgery of the vocal folds,” in Medical Imaging 2021: Image-Guided Procedures, Robotic Interventions, and Modeling, vol. 11598, p. 115981Z, 2021.
7. D. B. Camarillo, C. R. Carlson, J. K. Salisbury, Configuration tracking for continuum manipulators with coupled tendon drive. IEEE Trans. Robot. 25, 798–808 (2009).
8. M. P. Armacost, J. Adair, T. Munger, R. R. Viswanathan, F. M. Creighton, D. T. Crud, R. Sehra, Accurate and reproducible target navigation with the Stereotaxis Niobe® magnetic navigation system. J. Cardiovasc. Electrophysiol. 18, S26–S31 (2007).
9. A. K. Hoshiar, S. Jeon, K. Kim, S. Lee, J.-Y. Kim, H. Choi, Steering algorithm for a flexible microrobot to enhance guidewire control in a coronary angioplasty application. Micromachines 9, 617 (2018).
10. Johari, Shazlina & Shyan, L.. (2017). Stress-strain relationship of PDMS micropillar for force measurement application. EPJ Web of Conferences.
11. Shim, Sang & Yashin, Victor & Isayev, Avraam. (2004). Environmentally-friendly physico-chemical rapid ultrasonic recycling of fumed silica-filled poly(dimethyl siloxane) vulcanizate. Green Chemistry - GREEN CHEM. 
12. https://hypertextbook.com/facts/2004/MichaelShmukler.shtml
13. DAVIS, A. M. J., and K. B. RANGER. “A STOKES FLOW MODEL FOR THE DRAG ON A BLOOD CELL.” Quarterly of Applied Mathematics, vol. 45, no. 2, 1987, pp. 305–11
