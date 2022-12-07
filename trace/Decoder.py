import ctypes
import numpy as np
import os
from abc import ABCMeta, abstractmethod

class Decoder(metaclass=ABCMeta):
    def __init__(self,mod_order, num_sender, num_receiver, num_ofdm_sym) -> None:
        self.mod_order = mod_order
        self.num_sender = num_sender
        self.num_receiver = num_receiver
        self.num_ofdm_sym = num_ofdm_sym
    
    @abstractmethod
    def __call__(self,H,Y,w):
        pass
    
class ZFDecoder(Decoder):
    def __init__(self,mod_order, num_sender, num_receiver, num_ofdm_sym) -> None:
        super().__init__(mod_order, num_sender, num_receiver, num_ofdm_sym)
        
    def __call__(self,H,Y,w):
        X=np.zeros((self.num_ofdm_sym,self.num_sender),dtype=np.complex64)
        for i in range(self.num_ofdm_sym):
            X[i]=np.linalg.pinv(H[:,:,i]).T@Y[i,:]-w[:,i].T
        return X
            
    
    
        

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
    
# ND_POINTER_2D=np.ctypeslib.ndpointer(dtype=np.complex64,ndim=2,flags='C')
    

# path = os.path.dirname(__file__)
# cdll = ctypes.CDLL(os.path.join(path, "myDecoder.so"))
# my_decoder = cdll.my_decoder
# my_decoder.restype =ND_POINTER_2D
# my_decoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ND_POINTER_2D,ND_POINTER_2D]
