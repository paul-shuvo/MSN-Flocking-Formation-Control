#  MSN Flocking Formation Control

**This is an implementation of the MSN Flocking Formation Control presented in the paper [_Flocking for Multi-Agent Dynamic Systems: Algorithms and Theory_](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.121.7027&rep=rep1&type=pdf)**.

## Citation

### BibTex

```bib
@software{Shuvo_Kumar_Paul_MSN-Flocking-Formation-Control_An_implementation_2021,
author = {{Shuvo Kumar Paul}},
month = feb,
title = {{MSN-Flocking-Formation-Control: An implementation of multi-agent flocking formation control with specific formations that can follow a target without collision and can avoid obstacles}},
year = {2021}
}
```

### APA

```
Shuvo Kumar Paul. (2021). MSN-Flocking-Formation-Control: An implementation of multi-agent flocking formation control with specific formations that can follow a target without collision and can avoid obstacles [Computer software]
```

## Getting started

```shell
$ git clone
$ cd 
$ pip install
```
Run the python files for each of the five cases explained below.
```shell
$ python src/msn_1.py
```
## Project parameters:
- Number of sensor nodes: n =100.
- Space dimensions: m = 2.
- Desired distance among sensor node: d = 15.
- Scaling factor: k = 1.2 and interaction range r = k*d.
- Epsilon = 0.1 and Delta_t = 0.009 (These two parameters are optional and you can change them).

## Variable names
For easier understanding the variables are named similar to the equations described below.

For example,  <img src='alg_img\c_1_alpha.png' width=20px> is denoted as `C1_ALPHA` (all capital as it's a constant), <img src='alg_img\p_ik.png' width=30px> is denoted as `p_i_k`, and so on.

## Cases

### Case 1 - MSN Fragmentation

**Filename: msn_1.py**
- Randomly generates a connected network of 100 nodes in the area of 50x50. 
- Plots the initial deployment of the MSN of 100 nodes.
- Links the neighboringing nodes together by a line. 

#### Plots
- Plots the fragmentation of the nodes.
- Plots the velocity.
- Plots the connectivity.
- Plots the trajectory.

#### Algorithm 1:

<img src='alg_img\alg1.jpg' width=500px><br>
<img src='alg_img\alg11.jpg' width=500px><br>
<img src='alg_img\alg12.jpg' width=500px><br>
<img src='alg_img\alg13.jpg' width=500px><br>

#### Result
<img src='alg_img\msn_1.png' width=500px><br>


### Case 2 - Implements MSN Quasi-Lattice Formation with **static target**

**Filename: msn_2.py**

- Randomly generates a connected network of 100 nodes in the area of 50x50. 
- Sets up a target (gamma agent) as static point with its coordinate (x = 150, y =150). 
- Implements flocking behavior of the MSN.

#### Algorithm 2 for Case 2:

<img src='alg_img\alg2.jpg' width=500px><br>

#### Plots
- Plots the flocking of the nodes.
- Plots the velocity.
- Plots the connectivity.
- Plots the trajectory.

#### Result
<img src='alg_img\msn_2.png' width=500px><br>


### Case 3 - Implements MSN Quasi-Lattice Formation with **dynamic target (Sine wave trajectory)**

**Filename: msn_3.py**

- Randomly generates a connected network of 100 nodes in the area of 150x150. 
- Sets up a target (gamma agent) moving in a **sine wave trajectory.** 
- Implements flocking behavior of the MSN.

#### Algorithm 2 for Case 3:

<img src='alg_img\alg21.jpg' width=500px><br>

#### Plots
- Plots the flocking of the nodes following the target (Sine wave trajectory).
- Plots the velocity.
- Plots the connectivity.
- Plots the trajectory.
- Plots the center of mass and target trajectory.

#### Result
<img src='alg_img\msn_3.png' width=500px><br>

### Case 4 - Implements MSN Quasi-Lattice Formation with **dynamic target  (Circular trajectory)**

**Filename: msn_4.py**

- Randomly generates a connected network of 100 nodes in the area of 150x150. In this case you
- Sets up a target (gamma agent) moving in a **circular trajectory.** 
- Implements flocking behavior of the MSN.

#### Algorithm 2 for Case 3 (Same equation, change the values to make it a circular trajectory):

<img src='alg_img\alg21.jpg' width=500px><br>

#### Plots
- Plots the flocking of the nodes following the target (Circular trajectory).
- Plots the velocity.
- Plots the connectivity.
- Plots the trajectory.
- Plots the center of mass and target trajectory.

#### Result
<img src='alg_img\msn_4.png' width=500px><br>

### Case 5 - Implements MSN Quasi-Lattice Formation with **obstacle avoidance**

**Filename: msn_5.py**

- Randomly generates a connected network of 100 nodes in the area of 50x50. 
- Sets up a target (gamma agent) at the location of (200, 25). 
- Sets up an obstacle, circular in shape with a radius of 15 and its center location at (100,25). 
- Implements flocking behavior of the MSN avoiding obstacles.

#### Algorithm 3:

<img src='alg_img\alg3.jpg' width=500px><br>
<img src='alg_img\alg31.jpg' width=500px><br>
<img src='alg_img\alg32.jpg' width=500px><br>

#### Plots
- Plots the flocking of the nodes following the target, avoiding obstacles.
- Plots the velocity.
- Plots the connectivity.
- Plots the trajectory.
- Plots the center of mass and target trajectory.

#### Result
<img src='alg_img\msn_5.png' width=500px><br>
