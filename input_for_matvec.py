import h5py
import lattice_symmetries as ls
import numpy as np
import os
import yaml

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

def generate(basis_filename: str, output_filename: str, batch_size: int = 1):
    hamiltonian = load_hamiltonian(basis_filename)
    hamiltonian.basis.build()
    x = np.random.rand(hamiltonian.basis.number_states, batch_size) - 0.5
    x = np.asfortranarray(x)
    y = hamiltonian(x)

    with h5py.File(output_filename, "w") as out:
        out["/x"] = x.T
        out["/y"] = y.T


def main():
    generate("data/heisenberg_chain_10.yaml", "data/matvec/heisenberg_chain_10.h5", 1)
    generate("data/heisenberg_square_5x5.yaml", "data/matvec/heisenberg_square_5x5.h5", 1)

if __name__ == '__main__':
    main()
