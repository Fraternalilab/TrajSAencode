from TrajSAencode.SAlib import SADICT
from TrajSAencode.TrajCluster import TrajCluster

in_file = "/home/plasmids/testing/AllohubAnalysis/gamd_fbp_repl1_tc1.sasta"
cc = TrajCluster(SADICT, 4, None, 95)
cc.read_sasta(in_file)
cc.sa_traj = cc.sa_traj[::10]
#cc.read_sasta(in_file2)
cc.init_sim_mat()
cc.compute_similarity()
cc.cluster_traj()
cc.save_clusters("fbp_replO")
