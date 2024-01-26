import mdtraj as mdt
import argparse

def fix_pbc(traj):
    """Fixes the periodic boundary conditions of a trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory to fix.

    Returns
    -------
    traj : mdtraj.Trajectory
        Trajectory with fixed periodic boundary conditions.
    """
    anchors = [set(traj.topology.chain(i).atoms) for i in range(traj.n_chains)]

    new_traj = traj.image_molecules(
        inplace=False,
        anchor_molecules=anchors
    )

    return new_traj

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fixes the periodic boundary conditions of a trajectory."
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input trajectory file."
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output trajectory file."
    )

    args = parser.parse_args()

    traj = mdt.load(args.input)
    traj = fix_pbc(traj)
    traj.save_pdb(args.output)