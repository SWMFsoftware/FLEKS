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