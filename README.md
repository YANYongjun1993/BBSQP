# BBSQP
A MATLAB implementation of an integrated approach fot solving the mixed integer nonlinear programming (MINLP). This method is based on
Branch & Bound (BB) method, but does not require the relaxed nonlinear programming (NLP) at each node to be solved. Instead, early branching is
applied to reduce the computational burden of this numerical method. Then, a Sequential Quadratic Programming (SQP) is used to solve the relaxed NLP.
## Instructions
In the Test_case/ folder, two toy problems are used to calibrate the SQP and BB algorithm successively. For the SQP, the toy problem a spacecraft control problem with MPC framework.
Pls refer to the homework 4 of AS740 of UMICH for more derails. As for the BB, the first test problem in "Duran M A, Grossmann I E. An outer-approximation algorithm for a 
class of mixedinteger nonlinear programs[J]. Mathematical programming, 1986, 36(3): 307-339" is applied.

In the Eco-coasting folder, the proposed BBSQP method is used to solve a eco-coasting problem. Herein, fuel cut-off stragety is proposed to reveal the potential benefit of eco-coasting 
using the road grade preview. Pls refer to "Yana Y, Li N, Song Z, et al. Eco-Coasting Strategies Using Road Grade Preview: Evaluation and Online Implementation Based on Mixed Integer 
Model Predictive Control[J]. arXiv preprint arXiv:2111.07377, 2021." for more information about the problem formulation. 

Different computational acceleration strategies combination are applied to the solving process of the BBSQP solver for the eco-coasting problem. 
In this folder, branch on the best cost solution and the most fractional variable is chosen. In the Eco_coasting/Str1 folder, switch numbers constraints over the predictive horizon, infeasibility detection, and switch constraints are embedded in the branch and bound strategy to reduce the useless branch based on the physical knowledge of the plant.
In the Eco_coasting/Str2, warm start, infeasibility detection, and switch constraints are used.
## Contributing
FBstab is used as the QP solver. We must give credit to Dr. Dom's work to make this algorithm works[https://github.com/dliaomcp/fbstab-matlab].
