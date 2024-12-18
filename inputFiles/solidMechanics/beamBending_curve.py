import numpy as np


def convert_bulk_shear(bulk_modulus, shear_modulus):
    return 9.0 * bulk_modulus * shear_modulus / (3.0 * bulk_modulus + shear_modulus)


def max_bending(x, youngs_modulus, moment_inertia, beam_length, force):
    return (force * x * x) * (3.0 * beam_length - x) / (6.0 * youngs_modulus * moment_inertia)


def curve(**kwargs):
    x = np.squeeze(kwargs['totalDisplacement ReferencePosition trace'][0, :, 0])
    t = np.squeeze(kwargs['totalDisplacement Time'][:, 0])
    shear_modulus = 4.16667e9
    bulk_modulus = 5.5556e9
    beam_size = [80.0, 8.0, 4.0]
    traction = 1e6
    artificial_stiffness = 1.043

    # Calculate max deflection
    beam_I = beam_size[1] * (beam_size[2]**3) / 3
    youngs_modulus = convert_bulk_shear(bulk_modulus, shear_modulus)
    f = traction * beam_size[1] * beam_size[2]

    # In the example, the traction is scaled by time
    dy = np.zeros(np.shape(kwargs['totalDisplacement trace']))
    for ii, tb in enumerate(t):
        dy[ii, :, 1] = max_bending(x, youngs_modulus * artificial_stiffness, beam_I, beam_size[0], f * (tb + 1.0))

    return dy
