# Unit conversion

A physical quantity contains two parts: the value times the unit. For example, the weight of a table is $w = 50 kg = 50,000 g$, where 50 and 50,000 are the values, and kg and g are the units. If we define a new unit wg, which is defined as $1 wg = 2.5 kg$, then the table weight is $w = 20 wg$. In general, a physical quantity can be expressed as $u = \bar{u} u^*$, where $\bar{u}$ is the value/number and $u^*$ is the unit.  

## Mass normalization

With the CGS unit, the particle motion equation is
$$
\frac{dv}{dt} = \frac{q}{m}(E + \frac{v}{c} \times B) 
$$ 
In FLEKS, all quantities are normalized. The normalization should not introduce new constants, so we obtain $\frac{q_{cgs}^* E_{cgs}^* t_{cgs}^*}{m_{cgs}^* v_{cgs}^*} = 1$ from the equation above. 

From Gauss's law $\nabla\cdot E = 4\pi \rho$, we obtain $E_{cgs}^* = \rho_{cgs} ^*x_{cgs}^*$, where $\rho_{cgs}^* = q_{cgs}^*/x_{cgs}^{*3}$ is the charge density instead of mass density. 

From the expression above, $q_{cgs}^*=v_{cgs}^* \sqrt{m_{cgs}^* x_{cgs}^*}$ is obtained. 

In the code, $\bar{q}/(\bar{m}\bar{c}) = 1$ is required for a proton, so the normalization units $q_{cgs}^*$ and $m_{cgs}^*$ must satisfy $q_{cgs}^*/(m_{cgs}^*v_{cgs}^*) = q_{cgs,p}/(m_{cgs,p}c_{cgs})$, where $q_{cgs,p}$ and $m_{cgs,p}$ are proton charge and mass in CGS unit, respectively. Since $q_{cgs}^*=v_{cgs}^* \sqrt{m_{cgs}^* x_{cgs}^*}$, we obtain the normalization mass in CGS unit:
$$
m_{cgs}^* = x_{cgs}^*(\frac{c_{cgs}m_{cgs,p}}{q_{cgs,p}})^2.
$$ 
We note that $m_{cgs}^*$ also consists of the *value* and the unit *g*, and the expression above calculates its *value*. Its *value* in SI unit is calculated inside BATSRUS::ModPIC.f90, where the proton charge and mass are known in SI units. We need to convert the expression above into SI unit. Again, $m_{cgs}^*$ and $m_{SI}^*$ are the *values* of the normalization mass in CGS and SI unit, respectively. If the normalization mass is 1.5kg, then $m_{SI}^* = 1.5$ and $m_{SI}^* = 1500$. The conversion is $m_{cgs}^* = 1000m_{SI}^*$. For the charge unit, $1C = 3\times 10^9 esu$ assume the speed of light is $3 \times 10^8m/s$. So, we have the following expression:
$$
1000m_{SI}^* = 100*x_{SI}^*(\frac{100*c_{SI}*1000*m_{SI,p}}{3\times10^9 q_{SI,p}})^2 \\
m_{SI}^* = 10^7x_{SI}^*(\frac{m_{SI,p}}{q_{SI,p}})^2.
$$

Once $m_{SI}^*$ is obtained, it is passed to FLEKS and converted to $m_{cgs}^*$. 

## Length and velocity normalization

These two are free parameters. 

## Charge and current normalization

$q_{cgs}^*=v_{cgs}^* \sqrt{m_{cgs}^* x_{cgs}^*}$, as it is shown in previous section. 

$j_{cgs}^* = q_{cgs}^*v_{cgs}^*/x_{cgs}^{*3}$

## B and E normalization

B and E have the same unit in CGS. From previous section, we obtaied $\frac{q_{cgs}^* E_{cgs}^* t_{cgs}^*}{m_{cgs}^* v_{cgs}^*} = 1$. Substituting the expression of $q_{cgs}^*$ into it to obtain $E_{cgs}^* = B_{cgs}^* = \sqrt{\frac{m_{cgs}^*}{x_{cgs}^*3}}v_{cgs}^*$ 

# The consequence of changing light speed
Reducing the light speed changes the EM wave speed in the simulation, but it does change the interaction between magnetic field and particles. When we discuss the mass normalization, $\bar{q}/(\bar{m}\bar{c}) = 1$ instead of $\bar{q}/\bar{m} = 1$ is required. Because there is no assumption of $\bar{c}=1$. The light speed in the equation
$\frac{dv}{dt} = \frac{q}{m}(E + \frac{v}{c} \times B)$ is always the physical light speed no matter what $v_{cgs}^*$ is choosen. For example, if $v_{cgs}^* = 0.1c$, $\bar{c}$ has to be 10, otherwise, $\bar{q}/(\bar{m}\bar{c}) = 1$ can not be satisfied. 

Since reducing $v_{cgs}^*$ does not change the interaction between magnetic field and particle, the particle gyromotion, including gyro-radius and gyro-frequency, does not change. The inertial length is the gyro-radius of a particle with Alfven velocity, and neither gyromotion nor Alfven velocity changes, so the inertial length does not change with $v_{cgs}^*$.

How about the interaction between the electric field and particles? Since $\bar{c} = \bar{q}/\bar{m}$, which is not necessarily be 1, the Coulomb force for a proton should be $\bar{q}/\bar{m}E = \bar{c}E$. However, in the code, $\bar{c}$ is ignored, which suggests the Coulomb force in the simulation is $\bar{c}$ times weaker than in reality. Is this reasonable? What is its consequence? It is to be clarified. 

Most PIC simulations use reduced light speed, maybe the results are interpreted in a different way than FLEKS/MHD-EPIC, but I believe they also change the ratio between the Coulomb force and the $v \times B$ force. So, our approach of reducing light speed is probably reasonable. 

# Particle mover

## Boris algorithm

FLEKS uses Boris particle mover: 
$$
\frac{\bold{v}^{n+1} - \bold{v}^{n}}{\Delta t} = \frac{q}{m}(\bold{E}^{n+\theta} + \bar{\bold{v}}\times\bold{B})\\
\frac{\bold{x}^{n+1/2} - \bold{x}^{n-1/2}}{\Delta t} = \bold{v}^n
$$
where $\bold{v}$ is defined as $\bar{\bold{v}} = \frac{\bold{v}^{n+1} + \bold{v}^{n}}{2}$. The traditional method updates the velocity with the following steps.
* Step 1: acceleration 
    $$
    \bold{v}^- = \bold{v}^{n} + \frac{q\Delta t}{2m}\bold{E}^{n+\theta}
    $$
* Step 2: rotation
    $$ 
    \bold{a} = \bold{v}^- + \bold{v}^- \times \bold{t} \\
    \bold{v}^+ = \bold{v}^- + \bold{a} \times \bold{s}
    $$
    where $\bold{t} = \frac{q\Delta t}{2m} \bold{B}$, and $\bold{s} = \frac{2\bold{t}}{1+t^2}$.
* Step 3: acceleration
    $$
    \bold{v}^{n+1} = \bold{v}^+ + \frac{q\Delta t}{2m}\bold{E}^{n+\theta}
    $$

FLEKS solves the same velocity equation, but the implementation is slightly different.
* Step 1: acceleration 
    $$
    \bold{v}^- = \bold{v}^{n} + \frac{q\Delta t}{2m}\bold{E}^{n+\theta}
    $$
* Step 2: calculate $\bar{\bold{v}}$. With the expressions of $\bold{v}^+$ and $\bold{v}^-$ above, it is easy (it is probably not straightforwad, but I do not know how to number equations in markdown) to obtain $\bold{v}^{n+1} - \bold{v}^{n-1} = \bold{v}^+ - \bold{v}^- + \frac{q\Delta t}{m} \bold{E}^{n+\theta}$. Substitute it into the velocity update equation to remove electric field from the expression $\frac{\bold{v}^{+} - \bold{v}^{-}}{\Delta t} = \frac{q}{m}\bar{\bold{v}}\times\bold{B}$. Since $\bold{v}^+ = 2\bar{\bold{v}} - \bold{v}^-$, we obtain the equation for $\bar{\bold{v}}$:
  $$
  \bar{\bold{v}} = \bold{v}^- + \frac{q\Delta t}{2m} \bar{\bold{v}} \times \bold{B} = \bold{v}^- + \bar{\bold{v}} \times \bold{t}   
  $$
To solve this equation, we apply $\cdot \bold{t}$ to the expression of $\bar{\bold{v}}$ to obtain $\bar{\bold{v}} \cdot \bold{t} = \bold{v}^- \cdot \bold{t}$. Then we apply $\times \bold{t}$ to the expression of $\bar{\bold{v}}$ to obtain
$$
\bar{\bold{v}} \times \bold{t} = \bold{v}^- \times \bold{t} + \bar{\bold{v}} \times \bold{t} \times \bold{t} = \bold{v}^- \times \bold{t} + (\bar{\bold{v}} \cdot \bold{t})\bold{t} - t^2\bar{\bold{v}} = \bold{v}^- \times \bold{t} + (\bold{v}^- \cdot \bold{t})\bold{t} - t^2\bar{\bold{v}}
$$
Substitute the expression of $\bar{\bold{v}} \times \bold{t}$ to the expression of $\bar{\bold{v}}$ to obtain 
$$
\bar{\bold{v}} = \frac{\bold{v}^- + \bold{v}^- \times \bold{t} + (\bold{v}^- \cdot \bold{t})\bold{t}}{1+t^2}
$$
* Step 3: Then it is easy to obtain $\bold{v}^{n+1} = 2*\bar{\bold{v}} - \bold{v}^{n}$

## Relativistic Boris algorithm

With $\bold{u} = \gamma \bold{v}$, the velocity update equation is
$$
\frac{\bold{u}^{n+1} - \bold{u}^{n}}{\Delta t} = \frac{q}{m}(\bold{E}^{n+\theta} + \bar{\bold{u}}\times \frac{\bold{B}}{\gamma ^{n+1/2}})
$$
* Step 0: convert $\bold{v}$ to $\bold{u}$
* Step 1: acceleration
    $$
    \bold{u}^- = \bold{u}^{n} + \frac{q\Delta t}{2m}\bold{E}^{n+\theta}
    $$
* Step 2: calculate $\bar{\bold{u}}$ with
$$
\bar{\bold{u}} = \frac{\bold{u}^- + \bold{u}^- \times \bold{t} + (\bold{u}^- \cdot \bold{t})\bold{t}}{1+t^2}
$$
where $\bold{t} = \frac{q\Delta t}{2m} \bold{B}/\gamma^{n+1/2}$. From the definition of $\gamma$ and $\bold{u} = \gamma \bold{v}$, $\gamma ^{n+1/2} = \sqrt{1+(u^-/c)^2}$ is obtained. 
* Step 3: Then it is easy to obtain $\bold{u}^{n+1} = 2*\bar{\bold{u}} - \bold{u}^{n}$
* Step 4: Convert $\bold{u}$ back to $\bold{v}$