from dolfin import *; from mshr import *
import dolfin.common.plotting as fenicsplot
from matplotlib import pyplot as plt
# Define rectangular domain 
L = 4
H = 2

XMIN = 0.0; XMAX = L
YMIN = 0.0; YMAX = H

# Inflow velocity
VIN = 4.0

# Generate mesh (examples with and without a hole in the mesh) 
#resolution = 32
#mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), L*resolution, H*resolution)
#mesh = generate_mesh(Rectangle(Point(0.0,0.0), Point(L,H)) - Circle(Point(0.5,0.5*H),0.2), resolution)
#mesh = generate_mesh(Rectangle(Point(0.0,0.0), Point(L,H)) - Circle(Point(1.5,0.25*H),0.2) - Circle(Point(0.5,0.5*H),0.2) - Circle(Point(2.0,0.75*H),0.2), resolution)
resolution = 64
mesh = generate_mesh(Rectangle(Point(XMIN, YMIN), Point(L, H)) - Circle(Point(1.5, 0.25 * H), 0.2) 
                     - Circle(Point(0.5, 0.5 * H), 0.2) - Circle(Point(2.0, 0.75 * H), 0.2) 
                     - Circle(Point(3.0, 0.3 * H), 0.15) - Circle(Point(2.5, 0.6 * H), 0.1), resolution)


# Define periodic boundary condition
class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        # Upper boundary is "target domain" G
        return near(x[1], YMAX) and on_boundary

    # Map lower boundary (H) to upper boundary (G)
    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]+H    # Map lower boundary to upper boundary

# Create periodic boundary condition
pbc = PeriodicBoundary()

# Local mesh refinement (specified by a cell marker)
no_levels = 0
for i in range(0,no_levels):
    cell_marker = MeshFunction("bool", mesh, mesh.topology().dim())
    for cell in cells(mesh):
        cell_marker[cell] = False
        p = cell.midpoint()
        if p.distance(Point(0.5, 0.5 * H)) < 0.5:
            cell_marker[cell] = True
    mesh = refine(mesh, cell_marker)

plt.figure()
plot(mesh, title="Mesh")
#plt.savefig("Mesh.png")
plt.show()

# Generate mixed finite element spaces (for velocity and pressure)
VE = VectorElement("CG", mesh.ufl_cell(), 2)
QE = FiniteElement("CG", mesh.ufl_cell(), 1)
WE = VE * QE


W = FunctionSpace(mesh, WE, constrained_domain=pbc)
V = FunctionSpace(mesh, VE, constrained_domain=pbc)
Q = FunctionSpace(mesh, QE, constrained_domain=pbc)


# Define trial and test functions
w = Function(W)
(u, p) = (as_vector((w[0],w[1])), w[2])
(v, q) = TestFunctions(W) 


# Examples of inflow and outflow conditions

uin = Expression(("VIN*4*(x[1]*(YMAX-x[1]))/(YMAX*YMAX)", "0."), YMAX=YMAX, VIN=VIN, element = V.ufl_element()) 
#pout = 0.0

# Inflow boundary (ib), outflow boundary (ob) and wall boundary (wb)
ib = Expression("near(x[0],XMIN) ? 1. : 0.", XMIN=XMIN, element = Q.ufl_element())
ob = Expression("near(x[0],XMAX) ? 1. : 0.", XMAX=XMAX, element = Q.ufl_element()) 
#wb = Expression("x[0] > XMIN + DOLFIN_EPS && x[0] < XMAX - DOLFIN_EPS ? 1. : 0.", XMIN=XMIN, XMAX=XMAX, element = Q.ufl_element())
wb = Expression("x[0] > XMIN + DOLFIN_EPS && x[0] < XMAX - DOLFIN_EPS && x[1] > YMIN + DOLFIN_EPS && x[1] < YMAX - DOLFIN_EPS ? 1. : 0.", XMIN=XMIN, XMAX=XMAX, YMIN=YMIN, YMAX=YMAX, element = Q.ufl_element())


h = CellDiameter(mesh)
C = 1.0e10
gamma = C/h

f = Expression(("0.0","0.0"), element = V.ufl_element())

# Define variational problem on residual form: r(u,p;v,q) = 0
residual = ( - p*div(v)*dx + inner(grad(u), grad(v))*dx + div(u)*q*dx + 
            gamma*(ib*inner(u - uin, v) + wb*inner(u, v))*ds - inner(f, v)*dx )


# Solve algebraic system 
solve(residual == 0, w) 


#!rm results-NS/*

# volumetric flow rate Q
n = FacetNormal(mesh)
Q = assemble(dot(u, n) * ib * ds)
print("Volumetric Flow Rate (Q):", Q)

# Extract pressure values at inlet and outlet
p_values = w.sub(1).compute_vertex_values(mesh)  # Extract pressure values from the function w

# Find indices of vertices corresponding to inlet and outlet
inlet_vertex_indices = [vertex.index() for vertex in vertices(mesh) if near(vertex.point()[0], XMIN)]
outlet_vertex_indices = [vertex.index() for vertex in vertices(mesh) if near(vertex.point()[0], XMAX)]

# Compute pressure at inlet and outlet
p_inlet = p_values[inlet_vertex_indices[0]]  # Assuming there's only one vertex at the inlet
p_outlet = p_values[outlet_vertex_indices[0]]  # Assuming there's only one vertex at the outlet

# Calculate pressure drop
delta_p = p_outlet - p_inlet

# Calculate modulus of Q and modulus of delta P
Q_magnitude = sqrt(Q**2)
delta_p_magnitude = abs(delta_p)

# Compute the ratio of Q magnitude to delta P magnitude
ratio = Q_magnitude / delta_p_magnitude

# Output the result
print("Ratio of |Q| to |delta P|:", ratio)

# Define cross-sectional area A
A = H

# Calculate Darcy velocity q
q = Q / A

# Output Darcy velocity
print("Darcy Velocity (q):", q)
q1_magnitude = sqrt(q**2)

# Compute the ratio of q magnitude to delta P magnitude
q_ratio = q1_magnitude / delta_p_magnitude
print("Ratio of |q| to |delta P|:", q_ratio)

# Output pressure drop
print("(p_inlet):", p_inlet)
print("(p_outlet):", p_outlet)
print("Pressure Drop (Delta p):", delta_p)


# Plot solution
#plt.figure()
#plot(u, title="Velocity")
plt.savefig("velocity.png")
#plt.figure()
#plot(p, title="Pressure")
plt.savefig("pressure.png")        
#plt.show()

# Export files
#!tar -czvf results-Stokes.tar.gz results-NS
#files.download('results-Stokes.tar.gz')

# Plot pressure and ratio
plt.figure()
plt.plot(p_values, [ratio]*len(p_values), linestyle='-')
plt.xlabel('Pressure')
plt.ylabel('|Q| / |ΔP|')
plt.title('For V_in=4')
plt.grid(True)
plt.show()
