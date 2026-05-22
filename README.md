# Great-Crested-Grebe-Optimizer
A Seasonally Switched, Nature-Inspired Metaheuristic for Engineering Optimization
This repository provides a clean MATLAB implementation of the Great Crested Grebe Optimizer, abbreviated as GCGO, for continuous optimization problems.
GCGO is a seasonally switched nature-inspired metaheuristic optimizer. It models two plumage phases of the great crested grebe and uses four behavior-inspired mechanisms to coordinate global exploration and local exploitation.
## Main Mechanisms
The algorithm contains four core mechanisms:
1. Molting-based phase switching  
   The stochastic molting factor controls the transition between the winter plumage phase and the summer plumage phase.
2. Population density cycle  
   The active population size is dynamically adjusted during the iterative search process.
3. Winter plumage phase  
   Diving foraging and defense escape are used to enhance global exploration and population diversity.
4. Summer plumage phase  
   Courtship and reproduction and brood carrying are used to enhance coordinated exploitation and local refinement.
## File Structure
GCGO/
├── GCGO.m
├── main_demo.m
└── README.md
# Requirements
MATLAB R2020b or later
# Default Parameter Settings
| Parameter           |  Value | Description                        |
| ------------------- | -----: | ---------------------------------- |
| `molting_threshold` | 1.9290 | Molting phase switching threshold  |
| `eta_B`             |   0.03 | Population cycling amplitude ratio |
| `alpha_N`           |   0.05 | Logarithmic modulation coefficient |
| `beta`              |    1.5 | Levy flight parameter              |

# Citation
Gang Zheng, Jiange Kou, Zhiguo Yang, Xiangkai Shen, Yixuan Wang, and Yan Shi,
Great Crested Grebe Optimizer: A Seasonally Switched, Nature-Inspired Metaheuristic for Engineering Optimization.
