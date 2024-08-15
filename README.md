# add_fragment: A Python Module for Functionalizing a Host Molecule

Often, you may want to add a fragment to a host system, which could be a simple organic molecule, a cage structure, or a Metal-Organic Framework (MOF) or Covalent Organic Framework (COF). This need arises when you wish to study the effect of a functional group on the chemistry of the system or functionalize a periodic system at a specific site for further analysis. With the `add_fragment` module, you can now easily accomplish this task.

The logic behind the code is simple and straightforward. It identifies the site you wish to functionalize, then aligns and translates the fragment to that site. You can include as many sites as you need.

# Dependencies
The code requires Python >=3.9 but if you have a lower version edit the  pyproject.toml file with the version of Python you wish before installing. So far there are only two dependencies, ase and pytests so you can also use the version that matches your Python version. 

# Installation

Since the module is still in development, it is currently available only via Git installation. A PyPI installation will be provided at a later stage.

```
git clone git@github.com:bafgreat/add_fragment.git
cd add_fragment
pip install .
```

# USAGE

Using the code is straightforward. The following example shows how to functionalise host system.
```
from add_fragment.io import load_system
from add_fragment.insert_fragment import functionalise_host

host_system = load_system('filenames_of_host')
sub_fragment1 = load_system('filenames_of_fragment')
functionalised_system = functionalise_host(host_system, sub_fragment1)
functionalised_system.write('combined_system.xyz')

```
### NB
Note that both systems should have a dummy atom with symbol "X" to indicate the functionalisation site.
