#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke and Dr. Shabhajit Das"
__status__ = "production"
import numpy as np
from ase.neighborlist import NeighborList
from ase.data import covalent_radii


def find_third_point(p1, p2, bond_length):
    '''
    A function to find a third point. This is important
    to obtain a new point which you can use to rotate your
    fragment during alignment. The function takes two
    points and returns a new point in the direction of
    the first point.

    For instance if you know the position of P1 and P2 and
    will like to get the position of P3, this function
    can enable you to quickly obtain the coordinates of P3.
                R-P2-P1----P3
    parameters:
    -----------
    p1 : xyz coordinates in np.ndarray format
    p2 : xyz coordinates in np.ndarray format
    bond_length : parameter that determines the distance
    of point p3 from point p1
    returns:
    np.ndarray : coordinates of the new point P3
    '''
    d = p2 - p1
    unit_vector = d / np.linalg.norm(d)
    return p1 - bond_length * unit_vector


def kabsch(coordinates1, coordinates2):
    """
    A function that uses th Kabsch algorithm to find the
    rotation matrix that can be used to rotate and align
    a coordinates onto another.
    parameters:
    -----------
    coordinates1 : np.ndarray : coordinates of the first fragment
    coordinates2 : np.ndarray : coordinates of the second fragment
    returns:
    rotation matrix : np.ndarray : optimal rotation matrix that
    will align both coordinates
    """
    covariant_matrix = np.dot(np.transpose(coordinates1), coordinates2)
    left_singular_vectors, eigen_value, right_singular_vectors = np.linalg.svd(covariant_matrix)
    d = (np.linalg.det(left_singular_vectors) * np.linalg.det(right_singular_vectors)) < 0.0
    if d:
        eigen_value[-1] = -eigen_value[-1]
        left_singular_vectors[:, -1] = -left_singular_vectors[:, -1]
    rotation_matrix = np.dot(left_singular_vectors, right_singular_vectors)
    return rotation_matrix


def find_connected_atoms(atoms, atom_index):
    '''
    A function that searches for the indice of tha atom directly bonded to
    a given atom index. In our context, it returns the index of the atom
    directly bonded to the dummy atom X. This will enable us to obtain
    a second index from which will can use their coordinates to compute
    the third index.
    parameters:
    -----------
    atoms: ase atom
    atom_index: index of the atom for which we need to
    find the index of the connected atom
    returns:
    index of the connected atom
    '''
    cutoffs = covalent_radii[atoms.numbers] + 0.3
    nl = NeighborList(cutoffs=cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    indices, _ = nl.get_neighbors(atom_index)
    return indices


def align_fragment_orientation(host_system, sub_fragment, bond_length):
    '''
    A function that aligns the the fragment to be added on
    the the site in the host system to be functionalise.
    The function begins by searching the indices of the sites
    to be functionalised. In this case the site contain a dummy atom
    with symbol X.

    N.B. This function can be modified to use a different approach
    for searching functionalisation sites.

    Once the sites are identified, it the searches for the nearest
    neighbor so that we can have a points with at least two coordinates
    to perform the the alignment.

    In the host system we search for a 3rd point which will enable
    us to obtain optimal direction for alignment.

    parameters:
    -----------
    host_system : ase.Atoms : The ASE atoms object representing the
    host system.
    sub_fragment : ase.Atoms : The ASE atoms object representing the
    fragment to be added.
    bond_length : parameter that determines the distance of point p3
    from point p1
    returns:
    ase.Atoms : The ASE atoms object representing the aligned fragment.
    sub_x : index of dummy atom in  sub_fragment
    host_x : index of dummy atom in host_system
    sub_neigh: index of the atom in the fragment that is directly
    connected to the dummy
    '''
    host_x = [atom.index for atom in host_system if atom.symbol == 'X'][0]
    sub_x = [atom.index for atom in sub_fragment if atom.symbol == 'X'][0]

    sub_neigh = find_connected_atoms(sub_fragment, sub_x)[0]
    host_neigh = find_connected_atoms(host_system, host_x)[0]

    third_point = find_third_point(host_system.positions[host_x],
                                   host_system.positions[host_neigh], bond_length)

    host_x_coords = np.array([host_system.positions[host_x].tolist(),
                              third_point.tolist()])
    sub_x_coords = sub_fragment.positions[[sub_x, sub_neigh]]

    host_x_coords -= host_system[host_x].position
    sub_x_coords -= sub_fragment[sub_x].position

    rotation_matrix = kabsch(sub_x_coords, host_x_coords)
    sub_fragment.positions -= sub_fragment[sub_x].position
    sub_fragment.positions = np.dot(sub_fragment.positions, rotation_matrix)
    sub_fragment.positions += sub_fragment[sub_x].position

    return sub_fragment, sub_x, host_x,  sub_neigh


def translate_fragment(host_system, sub_fragment, sub_neigh, host_x):
    """
    A function that translates the aligned fragment such that the position of
    atom that is directly connected to its dummy atom should intersect with
    the position of the dummy atom in the host system. This ensures that a bond
    is formed bewteen the host system and the fragment.
    Parameters:
    -----------
    host_system : ase.Atoms : The ASE atoms object representing the
    host system
    sub_fragment (ase.Atoms): The ASE atoms object representing the fragment to
    be translated.
    sub_neigh (int): Index of the atom directly connected to the dummy
    atom in sub_fragment.
    host_x (int): Index of the dummy atom in host_system.
    Returns:
    sub_fragment (ase.Atoms): The ASE atoms object representing the
    aligned and translated fragment.
    """
    translation_vector = host_system[host_x].position - sub_fragment[sub_neigh].position
    sub_fragment.positions += translation_vector
    return sub_fragment


def functionalise_host(host_system, sub_fragment, bond_length=1.5):
    """
    A function that functionalises a fragment onto a host system by performing
    the necessary alignments and translations.
    The function begins by aligning the fragment to be added on the site in the
    host system to be functionalised.

    Once the fragment is aligned, it finds the nearest neighbor so that
    we can have a points with at least two coordinates to perform the the
    alignment.

    In the host system we search for a 3rd point which will enable us to obtain
    optimal direction for alignment.

    Finally, the function translates the aligned fragment such that the
    position of atom that is directly connected to its dummy atom should
    intersect with the position of the dummy atom in the host system. T
    his ensures that a bond is formed bewteen the host system and the fragment.

    The function then combines the functionalised fragment with the host system
    and returns the resulting system.

    Note: The function assumes dummy atoms are used to represent the sites to
    be functionalised.
    Parameters:
    -----------
    host_system (ase.Atoms): The ASE atoms object representing the host system.
    sub_fragment (ase.Atoms): The ASE atoms object representing the
    fragment to be aligned.
    bond_length (float): Desired bond length between the connected atoms.

    Returns:
    functionalised host system (ase.Atoms): New functionalised system.
    """
    aligned_sub_fragment, sub_x, host_x,  sub_neigh = align_fragment_orientation(host_system, sub_fragment, bond_length)
    translated_fragment = translate_fragment(host_system, aligned_sub_fragment, sub_neigh, host_x)
    sub_frag = [i.index for i in translated_fragment if i.index != sub_x]
    new_fragment = translated_fragment[sub_frag]

    host_frag = [i.index for i in host_system if i.index != host_x]
    host_system = host_system[host_frag]

    functionalised_host = host_system + new_fragment
    return functionalised_host
