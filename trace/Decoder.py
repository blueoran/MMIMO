import ctypes
import numpy as np
import os
class c_complex(ctypes.Structure):
    # Complex number, compatible with std::complex layout
    _fields_ = [("real", ctypes.c_double), ("imag", ctypes.c_double)]

    def __init__(self, pycomplex):
        # Init from Python complex
        self.real = pycomplex.real
        self.imag = pycomplex.imag

    def to_complex(self):
        # Convert to Python complex
        return self.real + (1.j) * self.imag
    
ND_POINTER_2D=np.ctypeslib.ndpointer(dtype=np.complex64,ndim=2,flags='C')
    

path = os.path.dirname(__file__)
cdll = ctypes.CDLL(os.path.join(path, "myDecoder.so"))
my_decoder = cdll.my_decoder
my_decoder.restype =ND_POINTER_2D
my_decoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ND_POINTER_2D,ND_POINTER_2D]
