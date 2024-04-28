def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_type', '--tao_type', type = str, default = 'bncg', help = 'TAO algorithm type')
	parser.add_argument('-tao_max_funcs', '--tao_max_funcs', type = int, default = 10000, help = 'TAO maximum functions evaluations')
	parser.add_argument('-tao_monitor', '--tao_monitor', action = 'store_true', help = 'TAO monitor')
	parser.add_argument('-ls', '--lagrange_s', type = float, default = 1.0, help = 'Lagrange multiplier for structural material')
	parser.add_argument('-lr', '--lagrange_r', type = float, default = 5.0, help = 'Lagrange multiplier for responsive material')
	parser.add_argument('-tao_ls_type', '--tao_ls_type', type = str, default = 'more-thuente', help = "TAO line search")
	parser.add_argument('-tao_view', '--tao_view', action = 'store_true', help = "View convergence details")
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 100, help = 'Number of TAO iterations')
	parser.add_argument('-tao_gatol', '--tao_gatol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is less than this')
	parser.add_argument('-tao_grtol', '--tao_grtol', type = float, default = 1.0e-7, help = 'Stop if relative norm of gradient is less than this')
	parser.add_argument('-tao_gttol', '--tao_gttol', type = float, default = 1.0e-7, help = 'Stop if norm of gradient is reduced by this factor')
	parser.add_argument('-vs', '--volume_s', type = float, default = 0.3, help = 'Volume percentage for structural material')
	parser.add_argument('-vr', '--volume_r', type = float, default = 0.3, help = 'Volume percentage for responsive material')
	parser.add_argument('-k', '--kappa', type = float, default = 1.0e-2, help = 'Weight of Modica-Mortola')
	parser.add_argument('-e', '--epsilon', type = float, default = 5.0e-3, help = 'Phase-field regularization parameter')
	parser.add_argument('-o', '--output', type = str, default = 'output1', help = 'Output folder')
	parser.add_argument('-m', '--mesh', type = str, default = 'motion.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-es', '--esmodulus', type = float, default = 0.01, help = 'Elastic Modulus for structural material')
	parser.add_argument('-er', '--ermodulus', type = float, default = 1.0, help = 'Elastic Modulus for responsive material')
	parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
	parser.add_argument('-q', '--power_q', type = float, default = 2.0, help = 'Power for multiple-well function')
	parser.add_argument('-s', '--steamy', type = float, default = 1.0, help = 'Initial stimulus')
	options = parser.parse_args()
	return options

options = parse()

from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)
VVVV = VectorFunctionSpace(mesh, 'CG', 1, dim = 4)

# Create initial design
rho =  Function(VVVV, name = "Design variable")
rho_i = Function(V, name = "Material density")
rhos = Function(V, name = "Structural material")  # Structural material 1(Blue)
rhor = Function(V, name = "Responsive material")  # Responsive material 2(Red)

sxx = Function(V, name = "Stimulus xx") # Right 
syy = Function(V, name = "Stimulus yy")	# Up

# Create initial design and stimulus
x, y = SpatialCoordinate(mesh)
mesh_coordinates = mesh.coordinates.dat.data[:]
M = len(mesh_coordinates)

rhos.interpolate(Constant(options.volume_s))
rhos.interpolate(Constant(1.0), mesh.measure_set("cell", 4))
rhor.interpolate(Constant(options.volume_r))
rhor.interpolate(Constant(0.0), mesh.measure_set("cell", 4))

sxx.interpolate(Constant(options.steamy))
syy.interpolate(Constant(options.steamy))

rho = as_vector([rhos, rhor, sxx, syy])
rho = interpolate(rho, VVVV)
###### End Initial Design #####

# Define the constant parameter used in the problem
kappa = Constant(options.kappa)
lagrange_r = Constant(options.lagrange_r)
lagrange_s = Constant(options.lagrange_s)

# Total volume of the domain |omega|
omega = assemble(interpolate(Constant(1.0), V) * dx)

delta = Constant(1.0e-6)
epsilon = Constant(options.epsilon)
kappa_d_e = Constant(kappa / epsilon)
kappa_m_e = Constant(kappa * epsilon)

# Define the predescribed displacements
fxx = Constant((-1.0, 0.0))
fyy = Constant((0.0, -1.0))

u_star_xx = Constant((1.0, 0.0)) # Right
u_star_yy = Constant((0.0, 1.0)) # Up

# Young's modulus of the beam and poisson ratio
E_v = Constant(delta)
E_s = Constant(options.esmodulus)
E_r = Constant(options.ermodulus)
nu = Constant(0.3) #nu poisson ratio

mu_v = E_v/(2 * (1 + nu))
lambda_v = (E_v * nu)/((1 + nu) * (1 - 2 * nu))

mu_s = E_s/(2 * (1 + nu))
lambda_s = (E_s * nu)/((1 + nu) * (1 - 2 * nu))

mu_r = E_r/(2 * (1 + nu))
lambda_r = (E_r * nu)/((1 + nu) * (1 - 2 * nu))

def v_v(rho):
	return 1 - rho.sub(0) - rho.sub(1)

def v_s(rho):
	return rho.sub(0)

def v_r(rho):
	return rho.sub(1)

# Define h_v(rho)=rho_v^(p)
def h_v(rho):
	return pow((1 - rho.sub(0) - rho.sub(1)), options.power_p)

# Define h_s(rho)=rho_s^(p)
def h_s(rho):
	return pow(rho.sub(0), options.power_p)

# Define h_r(rho)=rho_r^(p)
def h_r(rho):
	return pow(rho.sub(1), options.power_p)

def s_xx(rho):
	return rho.sub(2)

def s_yy(rho):
	return rho.sub(3)


# Define the double-well potential function
# W(x, y) = (x + y)^q * (1 - x)^q * (1 - y)^q
def W(rho):
	return pow((rho.sub(0) + rho.sub(1)), options.power_q) * pow((1 - rho.sub(0)), options.power_q) * pow((1 - rho.sub(1)), options.power_q)

# Define strain tensor epsilon(u)
def epsilon(u):
	return 0.5 * (grad(u) + grad(u).T)

# Define the residual stresses
def sigma_A(A, Id):
	return lambda_r * tr(A) * Id + 2 * mu_r * A

# Define the stress tensor sigma_v(u) for void
def sigma_v(u, Id):
	return lambda_v * div(u) * Id + 2 * mu_v * epsilon(u)

# Define the stress tensor sigma_s(u) for structural material
def sigma_s(u, Id):
	return lambda_s * div(u) * Id + 2 * mu_s * epsilon(u)

# Define the stress tensor sigma_r(u) for responsive material
def sigma_r(u, Id):
	return lambda_r * div(u) * Id + 2 * mu_r * epsilon(u)


# Update Lagrange multipliers

# Define test function and beam displacement
vxx = TestFunction(VV)
vyy = TestFunction(VV)

uxx = Function(VV, name = "xx-displacement")
uyy = Function(VV, name = "yy-displacement")

pxx = Function(VV, name = "Adjoint variable xx")
pyy = Function(VV, name = "Adjoint variable yy")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)

# Define the Modica-Mortola functional
func1 = kappa_d_e * W(rho) * dx

func2_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func2_sub2 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func2_sub3 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx

func2 = kappa_m_e * (func2_sub1 + func2_sub2 + func2_sub3)

P = func1 + func2

# Penalty for stimulus on void + structural material
func3_xx = pow(v_v(rho), 2) * pow(s_xx(rho), 2) * dx
func4_xx = pow(v_s(rhos), 2) * pow(s_xx(rho), 2) * dx

func3_yy = pow(v_v(rho), 2) * pow(s_yy(rho), 2) * dx
func4_yy = pow(v_s(rhos), 2) * pow(s_yy(rho), 2) * dx

S = func3_xx + func4_xx + func3_yy + func4_yy

# Objective function + Modica-Mortola functional + Penalty
Obj_xx = 0.5 * inner(uxx - u_star_xx, uxx - u_star_xx) * dx(4)
Obj_yy = 0.5 * inner(uyy - u_star_yy, uyy - u_star_yy) * dx(4)

J = Obj_xx + Obj_yy + P + S

# Volume fraction penalties
func5 = lagrange_s * v_s(rho) * dx
func6 = lagrange_r * v_r(rho) * dx

# Objective function + volume penalties
JJ = J + func5 + func6

# Define the weak form for forward PDExx
a_forward_v_xx = h_v(rho) * inner(1.0e3 * sigma_v(uxx, Id), epsilon(vxx)) * dx
a_forward_s_xx = h_s(rho) * inner(sigma_s(uxx, Id), epsilon(vxx)) * dx
a_forward_r_xx = h_r(rho) * inner(sigma_r(uxx, Id), epsilon(vxx)) * dx
a_forward_xx = a_forward_v_xx + a_forward_s_xx + a_forward_r_xx

# Define the weak form for forward PDEyy
a_forward_v_yy = h_v(rho) * inner(1.0e3 * sigma_v(uyy, Id), epsilon(vyy)) * dx
a_forward_s_yy = h_s(rho) * inner(sigma_s(uyy, Id), epsilon(vyy)) * dx
a_forward_r_yy = h_r(rho) * inner(sigma_r(uyy, Id), epsilon(vyy)) * dx
a_forward_yy = a_forward_v_yy + a_forward_s_yy + a_forward_r_yy

L_forward_xx = inner(fxx, vxx) * ds(8) + s_xx(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(vxx)) * dx
L_forward_yy = inner(fyy, vyy) * ds(8) + s_yy(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(vyy)) * dx

L_forward_xx_s = s_xx(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(vxx)) * dx
L_forward_yy_s = s_yy(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(vyy)) * dx

R_fwd_xx = a_forward_xx - L_forward_xx
R_fwd_yy = a_forward_yy - L_forward_yy

R_fwd_xx_s = a_forward_xx - L_forward_xx_s
R_fwd_yy_s = a_forward_yy - L_forward_yy_s

# Define the Lagrangian
a_lagrange_v_xx = h_v(rho) * inner(1.0e3 * sigma_v(uxx, Id), epsilon(pxx)) * dx
a_lagrange_s_xx = h_s(rho) * inner(sigma_s(uxx, Id), epsilon(pxx)) * dx
a_lagrange_r_xx = h_r(rho) * inner(sigma_r(uxx, Id), epsilon(pxx)) * dx
a_lagrange_xx = a_lagrange_v_xx + a_lagrange_s_xx + a_lagrange_r_xx

a_lagrange_v_yy = h_v(rho) * inner(1.0e3 * sigma_v(uyy, Id), epsilon(pyy)) * dx
a_lagrange_s_yy = h_s(rho) * inner(sigma_s(uyy, Id), epsilon(pyy)) * dx
a_lagrange_r_yy = h_r(rho) * inner(sigma_r(uyy, Id), epsilon(pyy)) * dx
a_lagrange_yy = a_lagrange_v_yy + a_lagrange_s_yy + a_lagrange_r_yy

L_lagrange_xx = inner(fxx, pxx) * ds(8) + s_xx(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(pxx)) * dx
L_lagrange_yy = inner(fyy, pyy) * ds(8) + s_yy(rho) * h_r(rho) * inner(sigma_A(Id, Id), epsilon(pyy)) * dx

R_lagrange_xx = a_lagrange_xx - L_lagrange_xx
R_lagrange_yy = a_lagrange_yy - L_lagrange_yy

R_lagrange = R_lagrange_xx + R_lagrange_yy
L = JJ - R_lagrange

# Define the weak form for adjoint PDE
a_adjoint_v_xx = h_v(rho) * inner(1.0e3 * sigma_v(vxx, Id), epsilon(pxx)) * dx
a_adjoint_s_xx = h_s(rho) * inner(sigma_s(vxx, Id), epsilon(pxx)) * dx
a_adjoint_r_xx = h_r(rho) * inner(sigma_r(vxx, Id), epsilon(pxx)) * dx
a_adjoint_xx = a_adjoint_v_xx + a_adjoint_s_xx + a_adjoint_r_xx

a_adjoint_v_yy = h_v(rho) * inner(1.0e3 * sigma_v(vyy, Id), epsilon(pyy)) * dx
a_adjoint_s_yy = h_s(rho) * inner(sigma_s(vyy, Id), epsilon(pyy)) * dx
a_adjoint_r_yy = h_r(rho) * inner(sigma_r(vyy, Id), epsilon(pyy)) * dx
a_adjoint_yy = a_adjoint_v_yy + a_adjoint_s_yy + a_adjoint_r_yy

L_adjoint_xx = inner(uxx - u_star_xx, vxx) * dx(4)
L_adjoint_yy = inner(uyy - u_star_yy, vyy) * dx(4)

R_adj_xx = a_adjoint_xx - L_adjoint_xx
R_adj_yy = a_adjoint_yy - L_adjoint_yy

# Beam .pvd file for saving designs
beam = File(options.output + '/beam.pvd')
dJdrhos = Function(V)
dJdrhor = Function(V)
rho_res = Function(V, name = "Responsive")
rho_str = Function(V, name = "Structural")
rho_void = Function(V, name = "Void")

dJdsxx = Function(V)
dJdsyy = Function(V)

stimulusxx = Function(V, name = "Stimulus xx")
stimulusyy = Function(V, name = "Stimulus yy")

N = M * 4
index_s = []
index_r = []
index_sxx = []
index_syy = []

for i in range(N):
	if (i%4) == 0:
		index_s.append(i)
	if (i%4) == 1:
		index_r.append(i)
	if (i%4) == 2:
		index_sxx.append(i)
	if (i%4) == 3:
		index_syy.append(i)

def FormObjectiveGradient(tao, x, G):

	# Print volume fraction of structural material
	volume_s = assemble(v_s(rho) * dx)/omega
	print("The volume fraction(Vs) is {}".format(volume_s))

	# Print volume fraction of responsive material
	volume_r = assemble(v_r(rho) * dx)/omega
	print("The volume fraction(Vr) is {}".format(volume_r))
	print(" ")

	i = tao.getIterationNumber()
	if (i%5) == 0:
		rho_i.interpolate(rho.sub(1) - rho.sub(0))

		stimulusxx.interpolate(rho.sub(2))

		stimulusyy.interpolate(rho.sub(3))

		rho_str.interpolate(rho.sub(0))
		rho_res.interpolate(rho.sub(1))
		rho_void.interpolate(1 - rho.sub(0) - rho.sub(1))
		solve(R_fwd_xx_s == 0, uxx, bcs = bcs)
		solve(R_fwd_yy_s == 0, uyy, bcs = bcs)
		beam.write(rho_i, stimulusxx, stimulusyy, rho_str, rho_res, rho_void, uxx, uyy, time = i)

	with rho.dat.vec as rho_vec:
		rho_vec.set(0.0)
		rho_vec.axpy(1.0, x)


	# Solve forward PDEs
	solve(R_fwd_xx == 0, uxx, bcs = bcs)
	solve(R_fwd_yy == 0, uyy, bcs = bcs)

	# Solve adjoint PDEs
	solve(R_adj_xx == 0, pxx, bcs = bcs)
	solve(R_adj_yy == 0, pyy, bcs = bcs)

	# Evaluate the objective function
	objective_value = assemble(J)
	print("The value of objective function is {}".format(objective_value))

	# Compute gradiet w.r.t rho2 and rho3 and s
	dJdrhos.interpolate(assemble(derivative(L, rho.sub(0))).riesz_representation(riesz_map="l2"))
	dJdrhos.interpolate(Constant(0.0), mesh.measure_set("cell", 4))

	dJdrhor.interpolate(assemble(derivative(L, rho.sub(1))).riesz_representation(riesz_map="l2"))
	dJdrhor.interpolate(Constant(0.0), mesh.measure_set("cell", 4))
	
	dJdsxx.interpolate(assemble(derivative(L, rho.sub(2))).riesz_representation(riesz_map="l2"))
	dJdsyy.interpolate(assemble(derivative(L, rho.sub(3))).riesz_representation(riesz_map="l2"))

	G.setValues(index_s, dJdrhos.vector().array())
	G.setValues(index_r, dJdrhor.vector().array())
	G.setValues(index_sxx, dJdsxx.vector().array())
	G.setValues(index_syy, dJdsyy.vector().array())

	f_val = assemble(L)
	return f_val

# Setting lower and upper bounds
# No need to enforce box constraints
lb = as_vector((0, 0, -1, -1))
ub = as_vector((1, 1, 1, 1))
lb = interpolate(lb, VVVV)
ub = interpolate(ub, VVVV)

with lb.dat.vec as lb_vec:
	rho_lb = lb_vec

with ub.dat.vec as ub_vec:
	rho_ub = ub_vec

# Setting TAO solver
tao = PETSc.TAO().create(PETSc.COMM_SELF)
tao.setType('bncg')
tao.setObjectiveGradient(FormObjectiveGradient, None)
tao.setVariableBounds(rho_lb, rho_ub)
tao.setFromOptions()

# Initial design guess
with rho.dat.vec as rho_vec:
	x = rho_vec.copy()

# Solve the optimization problem
tao.solve(x)
tao.destroy()

end = time.time()
print("\nExecution time (in seconds):", (end - start))
