import ctypes
import numpy as np
import os
from abc import ABCMeta, abstractmethod
from copy import deepcopy

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
        X=np.zeros((self.num_ofdm_sym,self.num_sender),dtype=np.complex128)
        for i in range(self.num_ofdm_sym):
            X[i]=np.linalg.pinv(H[:,:,i]).T@Y[i,:]-w[:,i].T
        return X
            
    
class CDecoder(Decoder):
    arr_1d=np.ctypeslib.ndpointer(dtype=np.double,ndim=1,flags='C_CONTIGUOUS')
    arr_2d=np.ctypeslib.ndpointer(dtype=np.double,ndim=2,flags='C_CONTIGUOUS')
    

    def __init__(self,mod_order, num_sender, num_receiver, num_ofdm_sym) -> None:
        super().__init__(mod_order, num_sender, num_receiver, num_ofdm_sym)
        os.system('make')
        self.lib=np.ctypeslib.load_library('mydecoder.so','./')
        self.lib.CDecoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,self.arr_1d,self.arr_1d,self.arr_1d]
        self.lib.CDecoder.restype=ctypes.POINTER(ctypes.c_double)
        
    def __call__(self,H,Y,w):
        
        X=np.zeros((self.num_ofdm_sym,self.num_sender),dtype=np.complex128)
        for i in range(self.num_ofdm_sym):
            x=self.lib.CDecoder(self.mod_order,self.num_sender,self.num_receiver,1,np.ascontiguousarray(H[:,:,i].reshape(-1)).view(np.double).astype(ctypes.c_double),np.ascontiguousarray(Y[i,:]).view(np.double).astype(ctypes.c_double),np.ascontiguousarray(w[:,i]).view(np.double).astype(ctypes.c_double))
            X[i,:]=np.array([x[j]+1j*x[j+1] for j in range(0,2*self.num_sender,2)])
        return X
''' 
class CDecoder(Decoder):
    class c_complex(ctypes.Structure):
        # Complex number, compatible with std::complex layout
        _fields_ = [("real", ctypes.c_double), ("imag", ctypes.c_double)]
        
        def __init__(self,comp):
            self.real=comp.real
            self.imag=comp.imag
        
            
        def cst_to_comp(self):
            return np.complex128(self.real+1j*self.imag)
            
        
        
        
    # arr_1d=np.ctypeslib.ndpointer(dtype=c_complex,ndim=1,flags='C_CONTIGUOUS')
    # arr_2d=np.ctypeslib.ndpointer(dtype=c_complex,ndim=2,flags='C_CONTIGUOUS')
    arr_1d=np.ctypeslib.ndpointer(dtype=np.double,ndim=1,flags='C_CONTIGUOUS')
    arr_2d=np.ctypeslib.ndpointer(dtype=np.double,ndim=2,flags='C_CONTIGUOUS')
    

    def __init__(self,mod_order, num_sender, num_receiver, num_ofdm_sym) -> None:
        super().__init__(mod_order, num_sender, num_receiver, num_ofdm_sym)
        os.system('make')
        self.lib=np.ctypeslib.load_library('mydecoder.so','./')
        
        # self.lib.CDecoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.POINTER(self.comp)),ctypes.POINTER(self.comp),ctypes.POINTER(self.comp)]
        # self.lib.CDecoder.restype=ctypes.POINTER(self.comp)
        self.lib.CDecoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,self.arr_1d,self.arr_1d,self.arr_1d]
        self.lib.CDecoder.restype=ctypes.POINTER(ctypes.c_double)
        
        # self.lib=ctypes.cdll.LoadLibrary(os.path.join(os.path.dirname(__file__),'mydecoder.so'))
        # self.lib.CDecoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,np.ctypeslib,ctypes.c_void_p,ctypes.c_void_p]
        # self.lib.CDecoder.restype=ctypes.c_void_p
        
    def __call__(self,H,Y,w):
        
        X=np.zeros((self.num_ofdm_sym,self.num_sender),dtype=np.complex128)
        for i in range(self.num_ofdm_sym):
            # Hi=np.array([[self.c_complex(h) for h in hh] for hh in H[:,:,i]])
            # Yi=np.array([self.c_complex(y) for y in Y[i,:]])
            # wi=np.array([self.c_complex(ww) for ww in w[:,i]])
            # print(Y[i,:])
            # print(np.ascontiguousarray(Y[i,:]).view(np.float64))
            # print(np.ascontiguousarray(Y[i,:].view(np.float32)))
            # print(np.ascontiguousarray(Y[i,:].view(np.float64)).astype(ctypes.c_float))
            # Hi=deepcopy(H[:,:,i])
            # Yi=deepcopy(Y[i,:])
            # wi=deepcopy(w[:,i])
            # X[i]=self.lib.CDecoder(self.mod_order,self.num_sender,self.num_receiver,1,Hi.ctypes.data_as(ctypes.POINTER(ctypes.POINTER(self.comp))),Yi.ctypes.data_as(ctypes.POINTER(self.comp)),wi.ctypes.data_as(ctypes.POINTER(self.comp)))
            # print(self.lib.CDecoder(self.mod_order,self.num_sender,self.num_receiver,1,np.ascontiguousarray(H),np.ascontiguousarray(Y[i,:]),np.ascontiguousarray(w[:,i])))
            # print(self.lib.CDecoder(self.mod_order,self.num_sender,self.num_receiver,1,np.ascontiguousarray(H[:,:,i]).view(np.float64).astype(ctypes.c_double),np.ascontiguousarray(Y[i,:]).view(np.float64).astype(ctypes.c_double),np.ascontiguousarray(w[:,i]).view(np.float64).astype(ctypes.c_double)))
            # print(Y[i,:])
            x=self.lib.CDecoder(self.mod_order,self.num_sender,self.num_receiver,1,np.ascontiguousarray(H[:,:,i].reshape(-1)).view(np.double).astype(ctypes.c_double),np.ascontiguousarray(Y[i,:]).view(np.double).astype(ctypes.c_double),np.ascontiguousarray(w[:,i]).view(np.double).astype(ctypes.c_double))
            X[i,:]=np.array([x[j]+1j*x[j+1] for j in range(0,2*self.num_sender,2)])
            # print(X[i,:])
            # for i in range(2*self.num_sender):
            #     print(x[i])
            # print(np.frombuffer(x,np.double,2*self.num_sender * np.dtype(np.double).itemsize))
        # print(X)
        # X=np.stack(X)
        return X
'''      


    
# ND_POINTER_2D=np.ctypeslib.ndpointer(dtype=np.complex128,ndim=2,flags='C')
    

# path = os.path.dirname(__file__)
# cdll = ctypes.CDLL(os.path.join(path, "myDecoder.so"))
# my_decoder = cdll.my_decoder
# my_decoder.restype =ND_POINTER_2D
# my_decoder.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ND_POINTER_2D,ND_POINTER_2D]
