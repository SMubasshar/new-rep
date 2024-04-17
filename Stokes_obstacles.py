#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from dolfin import *
import matplotlib.pyplot as plt

# Define rectangular domain
L = 4
H = 2

# Function to update obstacle positions based on current time `t`
def update_obstacle_positions(mesh, t):
    # Define new obstacle positions here based on time `t`
    # For example, move obstacles with a sinusoidal motion:
    x_shift = 0.1 * sin(t)  # Adjust the amplitude and frequency as needed
    obstacle_1_center = Point(1.5 + x_shift, 0.25 * H)
    obstacle_2_center = Point(0.5 + x_shift, 0.5 * H)
    obstacle_3_center = Point(2.0 + x_shift, 0.75 * H)
    
    # Create updated geometry for obstacles
    geometry = Rectangle(Point(0.0, 0.0), Point(L, H))              - Circle(obstacle_1_center, 0.2)              - Circle(obstacle_2_center, 0.2)              - Circle(obstacle_3_center, 0.2)
             
    # Generate new mesh
    mesh = generate_mesh(geometry, resolution)
    
    return mesh

# Generate initial mesh
resolution = 32
mesh = generate_mesh(Rectangle(Point(0.0, 0.0), Point(L, H)) - Circle(Point(1.5, 0.25 * H), 0.2) - Circle(Point(0.5, 0.5 * H), 0.2) - Circle(Point(2.0, 0.75 * H), 0.2), resolution)

# Create periodic boundary condition
class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool((near(x[0], 0) or near(x[0], L)) and on_boundary)

    def map(self, x, y):
        y[0] = x[0] - L if near(x[0], 0) else x[0] + L
        y[1] = x[1]

pbc = PeriodicBoundary()

# Generate function spaces
VE = VectorElement("CG", mesh.ufl_cell(), 2)
QE = FiniteElement("CG", mesh.ufl_cell(), 1)
WE = VE * QE

W = FunctionSpace(mesh, WE, constrained_domain=pbc)
V = FunctionSpace(mesh, VE, constrained_domain=pbc)
Q = FunctionSpace(mesh, QE, constrained_domain=pbc)

# Define trial and test functions
w = Function(W)
(u, p) = (as_vector((w[0], w[1])), w[2])
(v, q) = TestFunctions(W)

# Example of increased inflow pressure
uin = Expression(("4*(x[1]*(YMAX-x[1]))/(YMAX*YMAX) + 2", "0."), YMAX=H, element=V.ufl_element())

# Define variational problem on residual form
h = CellDiameter(mesh)
C = 1.0e10
gamma = C / h

f = Expression(("0.0", "0.0"), element=V.ufl_element())

residual = (-p * div(v) * dx + inner(grad(u), grad(v)) * dx + div(u) * q * dx +
            gamma * (ib * inner(u - uin, v) + wb * inner(u, v)) * ds - inner(f, v) * dx)

# Solve algebraic system
solve(residual == 0, w)

# Save solution to file
u1 = project(u, V)
p1 = project(p, Q)
file_u = File("results-Stokes/u.pvd")
file_p = File("results-Stokes/p.pvd")
file_u << u1
file_p << p1

# Plot solution
plt.figure()
plot(u1, title="Velocity")
plt.figure()
plot(p1, title="Pressure")

# Export updated mesh
updated_mesh = update_obstacle_positions(mesh, t=0.0)  # Pass the current time `t`
file_mesh = File("results-Stokes/mesh.pvd")
file_mesh << updated_mesh

# Plot updated mesh with obstacles
plot(updated_mesh, title="Mesh with obstacles", linewidth=0.5)
plt.show()

