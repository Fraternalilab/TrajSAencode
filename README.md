# TrajSAencode
Python version of the Structural alphabet encoding method

# Provisional Installation
The final release will be available as a package trough pip.

This package requires python 3.8 or newer

GNU scientific library 
### python packages
mdtraj

cffi

### Regular installation
Clone the package to an empty folder.
Go to the TrajSAencode folder and run:
```{py}
python setup.py install
```
### Local installation
In case one wants a local installation or does not have the permisions to install it with the other packages, one can run:
```{py}
python setup.py build_ext --inplace
```
Then add the path to the python path:
```{sh}
export PYTHONPATH="/location/of/TrajSAencode:$PYTHONPATH"
```
The pythonpath needs to be updated with every login so it is recommended to add it to the .bashrc

# Running the code
Encoding only pdbs
```{py}
python -m TrajSAencode.TrajEncode --pdb list_of_pdbs
```
This will create one sasta file (fasta like format) for every pdb.

Encoding trajectories
```{py}
python -m TrajSAencode.TrajEncode --pdb topology --traj list_of_trajectories
```
The topology can be both a pdb file or a .gro file.
The trajectories **must** match the topology.
The code will create one asta file for every trajectory.

The memory usage is low so one can use trajetories of any size. However, it is recommended to create an sliced trajectory without waters first to decrease the computational cost.


# Running the Clustering

A provisional example of how to run the clustering can be found on the script run_cluster.py
