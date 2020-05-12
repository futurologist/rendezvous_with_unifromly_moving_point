# rendezvous_with_unifromly_moving_point

Assume you have the initial positions and velocities ``xA0, vA0`` and ``xB0, vB0`` of spaceships A and B respectively. B moves with no acceleration and with constant velocity ``vB0``. Therefore, it travels uniformly along a straight line. Its motion is described as: ``xB = xB0 + t*vB0``. Spaceship A can turn on and off an acceleration of constant magnitude ``a0`` but can change its direction as it sees fit. The velocity of A should not exceed certain value ``v_max > 0``.

Since spaceship B travels uniformly, along a straight line with constant velocity ``vB0``, it actually defines an inertial coordinate system. In other words, if we translate the original coordinate system and attach it to B, the new system travels with constant velocity along a straight line and is therefore also inertial. The transformation is Galilean, so one can define the following change of coordinates (in both directions)

```matlab
y = x - xB0 - t*vB0
u = v - vB0

x = y + xB0 + t*vB0
v = u + vB0
```
in particular, for B for any moment of time ``t`` we get 

```matlab
yB = xB - xB0 - t*vB0 = xB0 + t*vB0 - xB0 - t*vB0  = 0``
```

At time ``t=0``,

```matlab
yA0 = xA0 - xB0  

uA0 = vA0 - vB0
```
We are trying to design the control in this new coordinate system and them move it back into the original one. So let us switch coordinates:
```matlab
y = x - xB
u = v - vB0
```

So in this new inertial coordinate system we are solving a problem of control theory and to engineer a good control, we would use as a Lyapunov function (a function that allows us to guarantee certain stable behavior and design the proper expression for the acceleration ``a``) the magnitude of the velocity squared ``L = norm(u)^2``. We want to design acceleration ``a`` so that the Lyapunov function in the initial phase of the motion monotonically and steadily decreases while the velocity reorients appropriately. 

Define the unit vector 
```matlab
L_unit = cross(x0A - xB0, v0A - vB0) / norm(cross(x0A - xB0, v0A - vB0))
``` 
Let in the coordinate system attached to B the motion of A satisfy the system of ordinary differential equations (these equations in both systems are Newton's, because both systems are inertial):
```matlab
dy/dt = u

du/dt = - a0 * (u - cross(L_unit, u)) / norm(u - cross(L_unit, u))
``` 
In other words, the acceleration is set to 
```matlab
a = - a0 * (u - cross(L_unit, u)) / norm(u - cross(L_unit, u))
```
Observe that by design ``norm(a) = a0``. Because the vectors ``u`` and ``cross(L_unit, u)`` are orthogonal and of equal magnitude (simply ``cross(L_unit, u)`` is the ninety degree rotation of vector ``u``), the denominator simplifies to 
```matlab
norm(u - cross(L_unit, u)) = sqrt( norm(u - cross(L_unit, u))^2 ) 
                           = sqrt( norm(u)^2 + norm(L_unit, u)^2 )
                           = sqrt( norm(u)^2 + norm(L_unit)^2*norm(u)^2 )
                           = sqrt( norm(u)^2 + norm(u)^2)
                           = sqrt(2) * norm(u)
```
So the system of differential equations simplifies to 
```matlab
dy/dt = u

du/dt = -(a0/sqrt(2)) * u/norm(u) + (a0/sqrt(2)) * cross(L_unit, u)) / norm(u)
``` 
The system is designed so that A always moves in the plane passing thorugh the origin B and perpendicular to the vector ``L_unit``. 

Because``u`` and ``cross(L_unit, u)`` are perpendicular, their dot product is 0, which allows us to calculate the time-derivative of the lyapunov function along the solutions to the system above (``u'`` means transpose of column-vector ``u``):
```matlab
d/dt( L ) = d/dt( norm(u)^2 ) = d/dt( u' * u ) = u' * du/dt 
          = u' * (  -(a0/sqrt(2)) * u/norm(u) 
            + (a0/sqrt(2)) * cross(L_unit, u)) / norm(u) )
          = -(a0/sqrt(2)) * u'*u / norm(u) 
            + (a0/sqrt(2)) * u'*cross(L_unit, u)) / norm(u) 
          = -(a0/sqrt(2)) * norm(u)^2 / norm(u)
          = -(a0/sqrt(2)) * norm(u)
          = - (a0/sqrt(2)) * sqrt(L)

d/dt( L ) = -(a0/sqrt(2)) * sqrt(L) < 0
```
which means that ``norm(u)`` decreases with time to 0, as desired.

The system of differential equations, that governs the motion, is non-linear and not explicitly solvable, so we have to integrate it numerically. To do that, I have chosen a geometric integrator method, where the system is split into two explicitly solvable systems, whose solutions are combined together to give (a very good approximation of) the solution to the original system. The systems are:
```matlab
dy/dt = u / 2

du/dt = -(a0/sqrt(2)) u / norm(u)
```  
and

```matlab
dy/dt = u / 2

du/dt = (a0/sqrt(2)) cross(L_unit, u) / norm(u)
```  
Initially, the second system is nonlinear, however after we calculate:
```matlab
d/dt(norm(u)*2) = d/dt (dot(u, u)) = 2 * dot(u, du/dt) 
                = 2 * dot(u, (a0/sqrt(2)) * cross(L_unit , u))
                = 2 * (a0/sqrt(2)) * dot(u, cross(L_unit , u)) 
                = 0
```
we conclude that during the motion defined by this system, the magnitude of the 
velocity is constant, i.e. ``norm(u) = norm(u0)`` where ``u0 = u(0)``. Thus, the systems, together with their solutions, now look like:
```matlab
First system:

dy/dt = u / 2

du/dt = -(a0/sqrt(2)) u / norm(u)

Solution:
y(t) = y0  +  h * u0/2  -  t^2 * a0 * u0 / (4*sqrt(2)*norm(u0));
u(t) = u - t * a0 * u0 / (sqrt(2)*norm(u0));
```  
and

```matlab
Second system:

dy/dt = u / 2

du/dt = (a0/(sqrt(2)*norm(u0))) cross(L_unit, u) 

Solution:
y(t) = y0 + (sqrt(2)*norm(u0)/a0) *( cross(L_unit, u0)
          + sin( t * a0/(sqrt(2)*norm(u0)) ) * u0  
          - cos( t  *a0/(sqrt(2)*norm(u0)) ) * cross(L_unit, u0) )     
u(t) = cos( t  *a0/(sqrt(2)*norm(u0)) ) * u0  
       + sin( t  *a0/(sqrt(2)*norm(u0)) ) * cross(L_unit, u0)
```  
The solution to the original system can be approximated as follows. Select a time step ``h``. Then if at time ``t`` the spaceship's position and velocity have been calculated to be ``y, u``, the updated spaceship's position and velocity at time ``t + h`` can be calculated by first letting the ship move along the solution of the second system starting from ``y, u`` for time ``h/2``, then move along the solution of the first system for time ``h`` and then move along the solution of the second system for time ``h/2``.   
