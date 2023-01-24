# cahn_hil
Solver of Cahn-Hilliard (CH) equation - an equation of mathematical physics which describes the process of phase separation.

# Requirements

Matlab R2018b or higher. 

# Description

The purpose of the Cahn-Hilliard theory is to provide an equation of motion of the composition field $u(\vec{r}, t)$. The normalized CH equation with the information about the current temperature $T$ can be written as:

```math
\frac{\partial u(\vec{r}, t)}{\partial t}=\nabla\left[m(u) \nabla\left(-\frac{T_{c}-T}{T_{c}}u+u^{3}-\nabla^{2} u\right)\right],
```

where $T_c$ is the critical temperature of phase separation and $\frac{T_{c}-T}{T_{c}}$ is the reduced temperature. For classical LLPS the mobility function $m(u)$ is constant. 

## Numerical simulation

The Cahn-Hilliard equation, in general, cannot be solved analytically. Thus, numerical simulation is required for the detailed comparison of the experiment with the theoretical behavior. In order to solve CH equation, it is discretized on a square lattice of $N_x \times N_y$ points. The boundary conditions are periodic. Considering lattice spacing as $h$, the space vector $\vec{r}$ can be rewritten as $\vec{r}=h \cdot(i, j)$, where $i \in[1, N_x]$ and $j \in[1, N_y]$. Thus, $u(\vec{r_{i,j}})\equiv u_{i, j}$ and the chemical potential: 
```math
\mu_{i, j}=\mu\left(u\left(\mathbf{r}_{i, j}\right)\right),
```
```math
\mu(u)=-\frac{T_{c}-T}{T_{c}}u+u^{3}-\nabla^{2} u.
```
The $\nabla m \nabla \mu$ can be rewritten via a second-order finite-difference formula:

```math
\nabla m \nabla \mu \sim \frac{1}{h^2}\left[m_{i+\frac{1}{2}, j}\left(\mu_{i+1, j}-\mu_{i, j}\right)+m_{i-\frac{1}{2}, j}\left(\mu_{i-1, j}-\mu_{i, j}\right)\right.\left.+m_{i, j+\frac{1}{2}}\left(\mu_{i, j+1}-\mu_{i, j}\right)+m_{i, j-\frac{1}{2}}\left(\mu_{i, j-1}-\mu_{i, j}\right)\right] .
```
The Laplacians, e.g. $\nabla^{2} u_{i,j}$ , which appear in $\mu(u)$ can be approximated by the five-point formula:

$$\nabla^{2} u_{i,j} \sim \frac{1}{h^{2}} \delta^{2}u_{i,j}= \frac{u_{i+1, j}+u_{i-1, j}+u_{i, j-1}+u_{i, j+1}-4 u_{i, j}}{h^{2}}.$$
The values of the mobility at the interstitial lattice points $\left(i \pm \frac{1}{2}, j\right)$ and \textcolor{white}{\-}$\left(i, j \pm \frac{1}{2}\right)$ can be approximated by linear interpolation:
$$
\begin{aligned}
m_{i \pm 1 / 2, j}=m\left[\frac{1}{2}\left(u_{i, j}+u_{i \pm 1, j}\right)\right], \textrm{\quad and \quad}  m_{i, j \pm 1 / 2}=m\left[\frac{1}{2}\left(u_{i, j}+u_{i, j \pm 1}\right)\right].
\end{aligned}
$$

After the discretization in space, the CH equations transforms into a system of $N_x \cdot N_y$ nonlinear differential equations:
```math
\frac{d \mathbf{u}(t)}{d t}=\mathbf{F}[\mathbf{u}(t)].
```
Here, $\mathbf{F}$ is the difference operator, representing the right-hand side of CH, and $\mathbf{u}(t)=\left(u_{1,1}(t),u_{2,1}(t), \ldots, u_{N_x, N_y}(t)\right)^{\top}$ is the vector of the normalized concentration values at the lattice points. The system can be integrated numerically using different schemes, provided in the current code. For example, the Euler solver:
```math
\mathbf{u}(t+\tau)=\mathbf{u}(t)+\frac{1}{2} \tau\mathbf{F}[\mathbf{u}(t)],
```
where $\tau$ is the time step. For the convergence of the solution, the stability of a consistent finite-difference scheme is required. To avoid the instability, we should observe the following inequality:
```math
\tau<\frac{h^{4}}{8 \gamma\left(8-\beta h^{2}\right)},
```
where $\gamma=m(u_0)$, $\beta=1-3 u_0^{2}$ and $u_0$ is the initial concentration.

The Cahn-Hilliard modeling is widely used for the description of phase separating systems such as alloys, or polymer blends and is an intense area of research. In addition, it is encountered in the literature for many other applications, such as planet formation and cancer growth. 

## Results

The current code simulates the evolution of the concentration 2D map with the possibility of creation .avi video and calculates the corresponding scattering pattern, which can be used for further analysis. The results are saved in the corresponding .mat files. The example of the simulated concentration map is presented below:

![u_0 3_0 1_0 006_15300_simple_resolution (1)](https://user-images.githubusercontent.com/82273226/214304391-c5737c6c-91a7-4e5b-a565-0f7f0cdedd69.gif)

The current code was used for the following articles:

* Reverse-engineering method for XPCS studies of non-equilibrium dynamics, **A. Ragulskaya**, V. Starostin, N. Begam, A. Girelli, H. Rahmann, M. Reiser, F. Westermeier, M. Sprung, F. Zhang, C. Gutt and F. Schreiber, IUCrJ 9, 439-448, 2022.

* Interplay between Kinetics and Dynamics of Liquid-Liquid Phase Separation in a Protein Solution Revealed by Coherent X-ray Spectroscopy, **A. Ragulskaya**, N. Begam, A. Girelli, H. Rahmann, M. Reiser, F. Westermeier, M. Sprung, F. Zhang, C. Gutt and F. Schreiber, Journal of Physical Chemistry Letters 12 (30), 7085-7090, 2021.



