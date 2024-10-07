# Build unit cell
from ase.build import bulk
cu = bulk("Cu", "fcc", a=3.6, cubic=True)
cu.set_calculator(calc)
cu1 = bulk("Cu", "fcc", a=3.6, cubic=True)
cu1.set_calculator(calc)

# Define presure function 
def get_pressure(cell):
    p = cell.get_stress(voigt=False)[0][0]+cu.get_stress(voigt=False)[1][1]+cu.get_stress(voigt=False)[2][2]
    return -p/3

# Define hydrostatic list
hydrostatic_strain = np.linspace(0.95,1.05,15)
# Get Energy, volume and pressure vector for those strains
E = np.zeros(len(hydrostatic_strain))
V = np.zeros(len(hydrostatic_strain))
P = np.zeros(len(hydrostatic_strain))
for i in range(len(hydrostatic_strain)):
    # Apply stress
    cell = cu.get_cell()
    cell *= hydrostatic_strain[i]
    cu1.set_cell(cell, scale_atoms=True) 
    # Find E[i] and P[i] 
    E[i] = cu1.get_potential_energy()/cu1.get_number_of_atoms()
    P[i] = get_pressure(cu1)
    V[i] = pow(np.array(cu1.cell)[0][0],3)/cu1.get_number_of_atoms()

# Plot pressure
plt.figure(0)
plt.plot(V,P,'r')
plt.xlabel("Volume [A^3]")
plt.ylabel("Pressure [GPa]")

# Those who browse the documentation of ASE can discover that it has functionality to fit the equation of state
# and extract parameters using more accurate methods. Compare your result above to those computed by ASE. 

# The call below takes two arrays: the variable V is an array of volumes/atom and E
# is an array of corresponding potential energy values, again per atom. If you want to execute this
# notebook cell, you will need to substitute your own arrays for V and E.
from ase.eos import EquationOfState
from ase.units import kJ

plt.figure(1)
eos = EquationOfState(V, E, eos="birchmurnaghan") # Birch-Murnaghan is a particular functional form fitted to the equation of state
v0, e0, B = eos.fit()
print('Bulk modulus: ', B / kJ * 1.0e24, 'GPa')
eos.plot()
annot_min(V,E)
