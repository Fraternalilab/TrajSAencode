from TrajSAencode.TrajLoader import TrajLoader
from TrajSAencode.SAEncoder import SAEncoder
from TrajSAencode.SAlib import SADICT
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
    arg = parser.parse_args()
    return arg


def main(args):
    pdbs = args.pdb
    trajectories = args.traj
    sa_encoder = SAEncoder(SADICT)
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
            sa_encoder.encode_protein(frame, name)
            nframes += 1
            print("Processed Frames: %s" % nframes, end="\r")

    sa_encoder.close_output()

    print("\nEncoding Finished")


if __name__ == "__main__":
    args = parse_args()
    main(args)

