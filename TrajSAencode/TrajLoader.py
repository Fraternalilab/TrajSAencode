# ===============================================================================
# Trajencode
# TrajLoader.py
# Class to hold and process the files to encode
# ===============================================================================

import mdtraj as md
import os


class TrajLoader:

    def __init__(self, topology, mdtrajectory, split_chains=False, chunk_size=1, start_f=0, skip=0, stride=None):
        """
        Holds the trajectory and extracts the C alphas and chains from it
        :param topology: name of the pdb to use as topology file
        :param mdtrajectory: name of the MD trajectory file
        :param split_chains: whether the chains should be considered different trajectories
        :param chunk_size: The number of frames to extract at the same time, 0 for all
        :param start_f: The starting frame number to put in the name of the output fasta string
        :param skip : number of frames to skip from the trajectory
        :param stride : Only read every stride-th frame
        """
        self.topology_file = topology
        self.mdtrajectory = mdtrajectory
        self.split_chains = split_chains
        self.chunk_size = chunk_size
        self.start_f = start_f
        self.skip = skip
        self.stride = stride
        self.topology, self.ca_indexes = self.read_topology()  # mdtraj object and list of CA indexes
        if split_chains:
            self.chains = self.get_chains()  # list of chains with their CA indexes
        else:
            self.chains = [self.ca_indexes]

    def read_topology(self):
        """
        Reads the topology and extracts the C alpha indexes
        :return: mdtraj trajectory and list of indexes
        """
        traj_top = md.load(self.topology_file)
        ca_indexes = traj_top.topology.select("name CA")
        return traj_top.atom_slice(ca_indexes), ca_indexes

    def get_chains(self):
        """
        Extracts the chains by looking at the residue sequence number
        :return: list of chains
        """
        chain = []
        old_resnum = -1
        chains = []

        for residue in self.topology.topology.residues:

            if residue.resSeq < old_resnum:
                chains.append(chain)
                chain = []
            chain.append(residue.atom(0).index)
            old_resnum = residue.resSeq

        chains.append(chain)
        return chains

    def frames(self):
        """
        Generator that yields the next frame and the name fasta name to put in the encoding
        :return:
        """
        # Check if topology file = md trajectory file. This means that the input is a pdb
        if self.topology_file == self.mdtrajectory:
            traj = self.topology
        else:
            traj = md.iterload(self.mdtrajectory, chunk=self.chunk_size, top=self.topology_file,
                               atom_indices=self.ca_indexes, skip=self.skip, stride=self.stride)

        # Iterate over the trajectory
        n_frames = self.start_f

        for chunk in traj:
            # Iterate over the number of chains
            for i, chain in enumerate(self.chains):
                sliced_traj = chunk.atom_slice(chain)
                sliced_frames = n_frames
                # Iterate over the frames for that chain
                for frame in sliced_traj:
                    sliced_frames += 1
                    yield frame.xyz[0], self.generate_name(i+1, sliced_frames)
            n_frames += self.chunk_size

    def generate_name(self, chain_num, frame_num):
        """
        Creates a fasta name for the encoding
        :param chain_num: The number of the chain
        :param frame_num: The number of the frame
        :return: Fasta string
        """
        top_name = os.path.split(self.topology_file)[1].split(".")[0]
        traj_name = os.path.split(self.mdtrajectory)[1].split(".")[0]
        return ">%s.%s.chain%s|%s" % (top_name, traj_name, chain_num, frame_num)











