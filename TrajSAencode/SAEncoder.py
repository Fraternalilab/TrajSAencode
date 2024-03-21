import numpy as np
from _encodeframe.lib import encode_frame
from cffi import FFI


class SAEncoder:

    def __init__(self, sa_dict):
        """
        Encodes trajectories using a library of fragments
        :param sa_dict: library of fragments to use
        """

        self.sa_dict = sa_dict
        self.fragment_size = int(len(self.sa_dict[list(self.sa_dict.keys())[0]]) / 3)
        # Convert the sa library to a format compatible with the C wrapper
        self.sa_library = np.ndarray(shape=(len(self.sa_dict) * self.fragment_size, 3),
                                     buffer=np.array([self.sa_dict[key] for key in sorted(self.sa_dict.keys())], dtype=np.float32),
                                     dtype=np.float32)
        self.sa_library /= 10
        # Generate a mapping of the SA fragment name to its index
        self.sa_code_map = self._generate_samap()
        self.output_file = {}
        self.output_file_name = ""
        self.ffi = FFI()


    def _generate_samap(self):
        """
        Generates a map of SA fragment names to their respective indexes
        :return:
        """
        samap = {}
        for i, key in enumerate(self.sa_dict.keys()):
            samap[i] = key
        return samap

    def _c_encode(self, frame):
        """
        Encodes the given frame, into a string of SA fragments
        :param frame: trajectory frame, (np.ndarray, shape = (number of residues,3))
        :return: string of the encoded frame
        """
        # Variables needed for the C function
        protein_length = int(frame.shape[0])
        n_windows = protein_length - self.fragment_size + 1
        n_fragments = int(len(self.sa_dict))
        encoding = np.zeros(n_windows, dtype=np.int32)
        mdframe = self.ffi.cast("float(*)[3]", frame.ctypes.data)
        fragments = self.ffi.cast("float(*)[3]", self.sa_library.ctypes.data)
        c_encoding = self.ffi.cast("int *", self.ffi.from_buffer(encoding, require_writable=True))
        # Call to the C function that does the encoding
        encode_frame(n_windows, n_fragments, self.fragment_size, mdframe, fragments, c_encoding)
        # Convert the list of indexes to SA fragment names
        econded_string = self._map_encoding(c_encoding, n_windows)
        return econded_string

    def _map_encoding(self, encoded_prot, n_windows):
        encoded_string = ""
        for i in range(n_windows):
            encoded_string += self.sa_code_map[encoded_prot[i]]
        return encoded_string

    def encode_protein(self, frame, name):
        """
        Handles open, closing and writing the output files
        :param frame: protein frame to encode
        :param name: fasta name of the frame to encode
        :return:
        """
        # If the file or chain has changed, close the file and open a new one
        if "%s.sasta" % name.split(">")[1].split("|")[0] != self.output_file_name:
            self.output_file_name = "%s.sasta" % name.split(">")[1].split("|")[0]
            if self.output_file_name not in self.output_file:
                self.output_file[self.output_file_name] = open(self.output_file_name, "w")
        # Encode the protein frame
        encoded_string = self._c_encode(frame)
        self.output_file[self.output_file_name].write(name)
        self.output_file[self.output_file_name].write("\n")
        self.output_file[self.output_file_name].write(encoded_string)
        self.output_file[self.output_file_name].write("\n")

    def close_output(self):
        for key in self.output_file:
            self.output_file[key].close()
        self.output_file = {}
