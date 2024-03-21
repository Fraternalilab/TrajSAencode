import numpy as np
import mdtraj as md

class TrajProcessor:

    def __init__(self, sa_dict, cutoff):
        """
        Extracts distances 
        """
        self.sa_dict = sa_dict
        self.cutoff = cutoff
        self.fragment_size = int(len(self.sa_dict[list(self.sa_dict.keys())[0]]) / 3)
        self.output_file = {}
        self.output_file_name = ""

    def compute_distances(self, frame, name):
        """
        Handles open, closing and writing the output files
        :param frame: protein frame to compute the distances
        :param name: fasta name of the frame to encode
        :return:
        """
        # If the file or chain has changed, close the file and open a new one
        if "%s_distances.out" % name.split(">")[1].split("|")[0] != self.output_file_name:
            self.output_file_name = "%s_distances.out" % name.split(">")[1].split("|")[0]
            if self.output_file_name not in self.output_file:
                self.output_file[self.output_file_name] = {} #open(self.output_file_name, "w")
        # Compute all distances from this frame
        for i,atom in enumerate(frame):
            ones = np.ones(len(frame))
            distances = np.linalg.norm(frame - atom, axis=1)
            ones[distances > self.cutoff] = 0.0
            if not i in self.output_file[self.output_file_name]:
                self.output_file[self.output_file_name][i] = ones
            else:
                self.output_file[self.output_file_name][i] += ones
        # also count frames
        if not "frames" in self.output_file[self.output_file_name]:
            self.output_file[self.output_file_name]["frames"] = 1
        else:
            self.output_file[self.output_file_name]["frames"] += 1

    def convert_to_fragments(self):
        """
        Turns distance count between CA's to counts between SA fragments
        Two fragments are considered in contact if at least one of its parts is in contact
        """ 
        for key in self.output_file:
            # key 1 is the file name and key 2 is the residues
            for key2 in self.output_file[key]:
                new_array = []
                if key2 == "frames":
                    continue
                # number of fragments we expect
                num_fragments = len(self.output_file[key][key2]) - (self.fragment_size - 1)
                if key2 >= num_fragments:
                    continue
                for i in range(0, num_fragments):
                    # number of times it was close enough
                    value = np.max(self.output_file[key][key2][i:i+self.fragment_size])
                    new_array.append(value)
                self.output_file[key][key2] = new_array


    def save_output(self):
        """
        Saves results to output file
        """
        for key in self.output_file:
            with open(key, "w") as out:
                for key2 in self.output_file[key]:
                    if key2 == "frames":
                        out.write("frames : %s\n" % self.output_file[key][key2])
                        continue
                    out.write("%s:[ " % key2)
                    for element in self.output_file[key][key2]:
                        out.write("%s " % element)
                    out.write("]\n")

