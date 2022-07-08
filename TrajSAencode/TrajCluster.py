import numpy as np
from _kabsch.lib import wrmsd_kabsch
from cffi import FFI
from sklearn.cluster import AgglomerativeClustering


class TrajCluster:
    def __init__(self, sa_dict, fragment_size, nclusters, max_dist):
        self.sa_dict = sa_dict
        self.sub_matrix = {}
        self.dis_matrix = []
        self.fragment_size = fragment_size
        self.ffi = FFI()
        self.create_sub_matrix()
        self.sa_traj = []
        self.clustering = None
        self.max_dist = max_dist
        self.nclusters = nclusters

    def create_sub_matrix(self):
        self.sub_matrix = {}
        keys = [key for key in sorted(self.sa_dict.keys())]
        for i in range(len(keys)):
            for j in range(i, len(keys)):
                if i == j:
                    self.sub_matrix.setdefault((keys[i], keys[j]), 0.0)
                else:
                    fragment1 = np.array(self.sa_dict[keys[i]], dtype=np.float32)
                    c_fragment1 = self.ffi.cast("float(*)[3]", fragment1.ctypes.data)
                    fragment2 = np.array(self.sa_dict[keys[j]], dtype=np.float32)
                    c_fragment2 = self.ffi.cast("float(*)[3]", fragment2.ctypes.data)
                    rmsd = wrmsd_kabsch(self.fragment_size, c_fragment1, c_fragment2)
                    self.sub_matrix.setdefault((keys[i], keys[j]), rmsd)
                    self.sub_matrix.setdefault((keys[j], keys[i]), rmsd)

    def read_sasta(self, file):
        with open(file, "r") as inn:
            for line in inn:
                line = line.rstrip()
                if not line.startswith(">"):
                    self.sa_traj.append(line)

    def clean_sasta(self):
        self.sa_traj = []

    def init_sim_mat(self):
        self.dis_matrix = np.zeros(shape=(len(self.sa_traj), len(self.sa_traj)))

    def compare_frames(self, frame1, frame2):
        dis = 0.0
        for k in range(len(frame1)):
            dis += self.sub_matrix[(frame1[k], frame2[k])]
        return dis

    def compute_similarity(self):
        for i, frame1 in enumerate(self.sa_traj):
            for j in range(i, len(self.sa_traj)):
                dis = self.compare_frames(frame1, self.sa_traj[j])
                self.dis_matrix[i][j] = dis
                self.dis_matrix[j][i] = dis
            print("Finished %s" % i)

    def cluster_traj(self):
        self.clustering = AgglomerativeClustering(n_clusters=self.nclusters,
                                                  distance_threshold=self.max_dist,
                                                  affinity="precomputed",
                                                  linkage="average").fit(self.dis_matrix)

    def save_clusters(self, name):
        files = {}
        nclusters = max(self.clustering.labels_) + 1
        for i in range(nclusters):
            files[i] = open("clust_%s_%s.sasta" % (name, i), "w")
        for i, frame in enumerate(self.sa_traj):
            clust = self.clustering.labels_[i]
            files[clust].write(">CLUST\n")
            files[clust].write("%s\n" % frame)
        for key in files:
            files[key].close()

    @staticmethod
    def read_clust_file(c_file):
        c_file_h = open(c_file, "r")
        mem = []
        for line in c_file_h:
            line = line.rstrip()
            if not line.startswith(">"):
                mem.append(line)
        c_file_h.close()
        return mem

    def expand_clusters(self, cluster_files):
        clust_dict = {}
        # Load cluster data
        for f in cluster_files:
            clust_dict[f] = {}
            clust_dict[f]["handler"] = open(f, "a")
            clust_dict[f]["data"] = self.read_clust_file(f)
        # Compare each frame with cluster data
        for i,frame in enumerate(self.sa_traj):
            min_clust = None
            min_val = 0.0
            for key in clust_dict:
                distances = []
                for frame2 in clust_dict[key]["data"]:
                    distances.append(self.compare_frames(frame, frame2))
                avg_dis = sum(distances)/len(distances)
                if (min_clust is not None and avg_dis < min_val) or min_clust is None:
                    min_val = avg_dis
                    min_clust = clust_dict[key]["handler"]
            # Add new element to cluster
            min_clust.write(">CLUST\n")
            min_clust.write("%s\n" % frame)
            print("Finished frame %s" % i)
        #close files
        for key in clust_dict:
            clust_dict[key]["handler"].close()















  

