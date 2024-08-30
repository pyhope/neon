import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
import matplotlib.pyplot as plt

# ******Important note******
# There is a bug in the MDAnalysis module where it cannot correctly handle periodic boundary conditions when the left box boundary (xlo, ylo, zlo) has a non-zero value (e.g., the first values of the BOX BOUNDS are 1.0897046722 for the .dump files in the directory). I converted the file to the LAMMPS data format using atomsk, which adjusted xlo, ylo, and zlo to 0 (and xhi, yhi, zhi, along with the coordinates of all atoms, were also modified accordingly), allowing the code to run correctly. Therefore, when running this code, please ensure that the left box boundary is set to 0.

# Define the parameters
max_num = 4  # the number of nearest neighbors to be considered
filename = '299000.lmp'  # input filename

# Load the lammps data file
u = mda.Universe(filename, format='DATA', atom_style='id type x y z') # dump format is also available: u = mda.Universe(filename, format='LAMMPSDUMP')
box = u.dimensions

# Shift the coordinates of neons to the center of the box (only for better visualization, does not affect the analysis)
u.atoms.positions += np.array([3.0, 3.0, 0.0])
u.atoms.wrap(box=box)

# Select the atoms of interest
Mg_and_O = u.select_atoms('type 1 or type 2') # Mg and O atoms
Ne = u.select_atoms('type 3') # Ne atoms

# Create a neighbor search tree
search_tree = AtomNeighborSearch(Mg_and_O, box=box) # box is for periodic boundary conditions

# Count the number of Mg and O atoms among the nearest neighbors of each Ne atom
Mg_counts = []
O_counts = []
vac_type = []

for atom in Ne:
    neighbors = search_tree.search(atom, 5.0)
    # neighbors = neighbors.subtract(atom)
    sorted_neighbors = sorted(neighbors, key=lambda x: np.linalg.norm(x.position - atom.position))[:max_num]
    # print(atom.index + 1)
    # for neighbor in sorted_neighbors:
    #     print(neighbor.index, neighbor.type, np.linalg.norm(neighbor.position - atom.position))

    Mg_count, O_count = 0, 0
    for neighbor in sorted_neighbors:
        if neighbor.type == '1':
            Mg_count += 1
        elif neighbor.type == '2':
            O_count += 1

    Mg_counts.append(Mg_count)
    O_counts.append(O_count)
    if Mg_count == 0 and O_count == max_num:
        vac_type.append(1)
    elif Mg_count == max_num and O_count == 0:
        vac_type.append(2)
    else:
        vac_type.append(1.5)

for i, atom in enumerate(Ne):
    print(f"Ne #{atom.index + 1}:")
    print(f"{Mg_counts[i]} Mg atoms, {O_counts[i]} O atoms among the nearest {max_num} neighbors")
    if Mg_counts[i] == 0:
        print("This Ne atom is occupying Mg vacancy")
    elif O_counts[i] == 0:
        print("This Ne atom is occupying O vacancy")
    else:
        print("This Ne atom is occupying interstitial site")

# Visualization
pos = Mg_and_O.positions
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=Mg_and_O.types.astype(int), lw=0, cmap='brg', vmin=1, vmax=2, s=160, depthshade=True)
pos = Ne.positions
ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=vac_type, cmap='brg', vmin=1, vmax=2, lw=3, edgecolor='violet', s=400, depthshade=False)
print(vac_type)
ax.set_xlabel('x (Å)')
ax.set_ylabel('y (Å)')
ax.set_zlabel('z (Å)')
ax.view_init(elev=14, azim=-32)
plt.show()
