import numpy as np


def mat2vec(R):
    """
    Convert a matrix to a 1D array by flattening it row-wise.
    """
    return np.asarray(R).flatten(order='C')


def vec2mat(vect, larg):
    """
    Create a matrix from a 1D array, arranging its elements row-wise.
    Completes the last row with zeros if necessary.
    """
    vect = np.asarray(vect)
    if larg <= 0:
        raise ValueError("larg must be a positive integer.")

    nv = len(vect)
    n = int(np.ceil(nv / larg))

    mat = np.zeros((n, larg), dtype=vect.dtype)
    mat_flat = mat.ravel(order='C')
    mat_flat[:nv] = vect
    return mat


def bytes2bits(s, dtype='uint8'):
    """
    Convert a list of numbers or a single number to a bit array.
    """
    if np.isscalar(s):
        s = [s]

    s = np.asarray(s, dtype=dtype)
    s_bytes = s.view(np.uint8)

    s_bits = np.unpackbits(s_bytes, bitorder='little')
    return s_bits.astype(np.uint8)


def bits2bytes(xbits, dtype='uint8'):
    """
    Convert a bit array to bytes.
    """
    xbits = np.asarray(xbits).flatten()

    if np.any((xbits != 0) & (xbits != 1)):
        raise ValueError("xbits must contain only 0 or 1.")

    if len(xbits) % 8 != 0:
        raise ValueError("Number of bits must be a multiple of 8.")

    bytes_array = np.packbits(xbits.astype(np.uint8), bitorder='little')

    target_dtype = np.dtype(dtype)
    if bytes_array.nbytes % target_dtype.itemsize != 0:
        raise ValueError(
            f"Packed byte size ({bytes_array.nbytes}) is not compatible with dtype={dtype}."
        )

    return bytes_array.view(target_dtype)