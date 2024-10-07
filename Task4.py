# Shear modulus
# Build unit cell
import ase
from ase.units import GPa
from ase.build import bulk
cu_s = bulk("Cu", "fcc", a=3.6, cubic=True)
cu_s.set_calculator(calc)

# Apply shear stress function 

cu_m = np.array(cu_s.cell.copy())
cu_m[0][1]=0.01*3.6
cu_s.set_cell(cu_m,scale_atoms=True)

s = cu_s.get_stress(voigt=False)[0][1]/0.01
print("Shear modulus:",s/GPa,'GPa')


# Poisson ratio
from ase.build import bulk
cu_p = bulk("Cu", "fcc", a=3.6, cubic=True)
cu_p.set_calculator(calc)

# e_x strain value for x directions
e_x = - 0.01
e_yzs = np.linspace(0,0.01,100)
y_stress = []

for e_yz in e_yzs:
    # modify cell
    cu_mp = np.array(cu_p.cell)
    cu_mp[0][0]=3.6*(1+e_x)
    cu_mp[1][1]=3.6*(1+e_yz)
    cu_mp[2][2]=3.6*(1+e_yz)
    cu_p.set_cell(ase.cell.Cell(cu_mp),scale_atoms=True)
    stress = cu_p.get_stress(voigt=False)
    y_stress.append(stress[1][1])

plt.plot(e_yzs,y_stress)
plt.plot(e_yzs,np.zeros(len(e_yzs)),'--')

y_stress_abs = np.abs(y_stress)
i = np.argmin(y_stress_abs)

print('Poisson ratio:',-e_yzs[i]/e_x)
