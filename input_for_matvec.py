import h5py
import lattice_symmetries as ls
import numpy as np
import os
import time
import yaml

np.random.seed(42)

def load_hamiltonian(filename: str):
    _, extension = os.path.splitext(filename)
    if extension != ".yaml" and extension != ".yml":
        raise ValueError(
            "Do not know how to read the basis from files with extension '{}'..."
            " Try specifying a path to a YAML file instead.".format(extension)
        )

    with open(filename, "r") as f:
        config = yaml.load(f, Loader=yaml.SafeLoader)
    basis = ls.SpinBasis.load_from_yaml(config["basis"])
    return ls.Operator.load_from_yaml(config["hamiltonian"], basis)


def generate(basis_filename: str, output_filename: str, batch_size: int = 1, skip: bool = False):
    hamiltonian = load_hamiltonian(basis_filename)
    tick = time.time()
    hamiltonian.basis.build()
    tock = time.time()
    print("{} spent building the basis using OpenMP".format(tock - tick))
    x = np.random.rand(hamiltonian.basis.number_states, batch_size) - 0.5
    if skip:
        return

    x = np.asfortranarray(x)
    tick = time.time()
    y = hamiltonian(x)
    tock = time.time()
    print("{} spent in matrix-vector using OpenMP".format(tock - tick))

    with h5py.File(output_filename, "w") as out:
        out["/representatives"] = hamiltonian.basis.states
        out["/x"] = x.T
        out["/y"] = y.T


def main():
    # generate("data/heisenberg_chain_10.yaml", "data/matvec/heisenberg_chain_10.h5", 1)
    # generate("data/heisenberg_square_5x5.yaml", "data/matvec/heisenberg_square_5x5.h5", 1)
    # generate("data/heisenberg_square_6x6.yaml", "data/matvec/heisenberg_square_6x6.h5", 1)
    # generate("data/old/heisenberg_chain_20.yaml", "data/matvec/heisenberg_chain_20.h5", 1)
    for i in [10, 12, 16, 20, 24, 28, 32]:
        generate(
            "data/old/heisenberg_chain_{}.yaml".format(i),
            "data/matvec/heisenberg_chain_{}.h5".format(i),
            1,
            skip=True
        )
    generate(
        "data/old/heisenberg_chain_24_symm.yaml",
        "data/matvec/heisenberg_chain_24_symm.h5",
    )
    generate(
        "data/old/heisenberg_chain_32_symm.yaml",
        "data/matvec/heisenberg_chain_32_symm.h5",
    )
    generate(
        "data/old/heisenberg_chain_36_symm.yaml",
        "data/large-scale/matvec/heisenberg_chain_36_symm.h5",
    )
    generate(
        "data/old/heisenberg_kagome_12_symm.yaml",
        "data/matvec/heisenberg_kagome_12_symm.h5",
    )


if __name__ == "__main__":
    main()
