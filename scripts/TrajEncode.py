from TrajSAencode.TrajLoader import TrajLoader
from TrajSAencode.SAEncoder import SAEncoder
from TrajSAencode.SAlib import SADICT
from TrajSAencode.TrajProcessor import TrajProcessor
import argparse


class InputError(Exception):
    """Errors in the input provided"""
    pass


def parse_args():
    parser = argparse.ArgumentParser(description="Encodes Trajectories and PDBs into Structural Alphabet (SA) strings,"
                                                 "multiple files can be provided at the same time using wildcards")
    parser.add_argument('--pdb', type=str, nargs='+', required=True, help="Input PDBs")
    parser.add_argument('--traj', type=str, nargs='+', required=False, default="", help="Input trajectories")
    parser.add_argument('--split', type=bool, required=False, default=True, help="whether or not to split chains")
    parser.add_argument('--chunk', type=int, required=False, default=1, help="Number of frames to load at once into memory")
    parser.add_argument('--start', type=int, required=False, default=0, help="Starting number for the output numbering")
    parser.add_argument('--skip', type=int, required=False, default=0, help="Frames to skip")
    parser.add_argument('--stride', type=int, required=False, default=1, help="stride the trajectory")
    parser.add_argument('--mode', type=str, required=False, default="all", help="Way to process the trajectory: Options: encode, distance, all")
    parser.add_argument('--cutoff', type=float, required=False, default=10.0, help="cutoff to use to pick which atoms to account when computing distances")
    arg = parser.parse_args()
    return arg


def main(args):
    pdbs = args.pdb
    cutoff = args.cutoff / 10
    trajectories = args.traj
    if args.mode in ["all", "encode"]:
        sa_encoder = SAEncoder(SADICT)
    if args.mode in ["all", "distance"]:
        traj_pros = TrajProcessor(SADICT, cutoff)
    # Some checks
    if len(pdbs) == 0:
        raise InputError("PDB file not found")
    elif len(pdbs) == 1:
        if len(trajectories) == 0:
            trajectories = pdbs
        pdb = pdbs[0]
    elif len(pdbs) > 1:
        if len(trajectories) > 0:
            raise InputError("When processing trajectory files only one pdb can be provided")
        else:
            trajectories = pdbs
            pdb = pdbs[0]
    if not pdb.split(".")[-1] in ["pdb", "gro"]:
        raise InputError("PDB files must end with .pdb or .gro")

    if len(pdbs) == 1:
        print("Topology extracted from: %s" % pdbs[0])
    print("Processing files: %s" % ", ".join(trajectories))

    # Iterate over the trajectories to encode
    for traj in trajectories:
        print("Processing file: %s\n" % traj)
        # If all are pdbs the trajectory and the topology are the same file
        if len(pdbs) > 1:
            pdb = traj
        # Process trajectory frame by frame
        traj_loader = TrajLoader(pdb, traj, split_chains=args.split, chunk_size=args.chunk, start_f=args.start,
                                 skip=args.skip, stride=args.stride)

        nframes = 0
        for frame, name in traj_loader.frames():
            if args.mode in ["all", "encode"]:
                sa_encoder.encode_protein(frame, name)
            if args.mode in ["all", "distance"]:
                traj_pros.compute_distances(frame, name)
            else:
                print("Mode not found")
                exit(1)
            nframes += args.chunk
            print("Processed Frames: %s" % nframes, end="\r")

    if args.mode in ["all", "encode"]:
        sa_encoder.close_output()
    if args.mode in ["all", "distance"]:
        traj_pros.convert_to_fragments()
        traj_pros.save_output()

    print("\nEncoding Finished")


if __name__ == "__main__":
    args = parse_args()
    main(args)

