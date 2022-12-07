
import numpy as np
from channel_analysis import *
import h5py
import commpy as cpy
from Decoder import *


def generate_samples(from_file,num_sample):
    h5log=h5py.File(from_file,'r')
    f=h5py.File(f'./sample_traces/sample_{num_sample}.hdf5',"w")
    f.attrs['samples_per_user']=h5log.attrs['samples_per_user']
    f.attrs['num_mob_ant']=h5log.attrs['num_mob_ant']
    f.attrs['frame_length']= h5log.attrs['frame_length']
    if 'measured_noise' in h5log.attrs.keys():
        f.attrs['measured_noise']= h5log.attrs['measured_noise']
    f['Pilot_Samples']=h5log['Pilot_Samples'][:num_sample,:,]
    f.close()
    h5log.close()
    



def calc_martix(userCSI,user,frame,method='zf'):
    user = 0
    frame=10
    lts_iq = h5log['Pilot_Samples'][frame,0,user*samps_per_user:(user+1)*samps_per_user,0]*1.+h5log['Pilot_Samples'][frame,0,user*samps_per_user:(user+1)*samps_per_user,1]*1j
    lts_iq /= 2**15
    offset = lts.findLTS(lts_iq)[0][0]+32
    print("LTS offset for user %d, frame %d: %d" % (user, frame, offset))
    if method=='conj':
        bws = np.transpose(np.conj(userCSI[frame,:,:,:]),(1,0,2))
        return bws	
    elif method=='zf':	
        bws = np.empty((userCSI.shape[2],userCSI.shape[1],userCSI.shape[3]),dtype='complex128')	
        for sc in range(userCSI.shape[3]):
            bws[:,:,sc] = np.linalg.pinv(userCSI[frame,:,:,sc])
        return bws
    else:
        raise ValueError('method must be conj or zf')



def simulate(decoder_class,mod_num,bws,noise,gen_len,binom_prob=0.5,apply_noise=False):
    if mod_num<=8:
        Modulation_type=cpy.PSKModem(mod_num)
    else:
        Modulation_type=cpy.QAMModem(mod_num)
    decoder=decoder_class(mod_num,bws.shape[0],bws.shape[1],bws.shape[2])
    msg_send = np.random.binomial(n=1,p=binom_prob,size=(bws.shape[2],bws.shape[0],int(np.log2(mod_num)*(1+gen_len//np.log2(mod_num)))))
    msg_send_mod=np.array([[Modulation_type.modulate(msg_send[i,j,:]) for j in range(msg_send.shape[1])] for i in range(msg_send.shape[0])])
    # print(np.sum(msg_send!=Modulation_type.demodulate(msg_send_mod.reshape(-1), 'hard').reshape(msg_send.shape)))
    noise=noise.astype('complex128')
    noise_random=np.zeros_like(msg_send_mod,dtype='complex128')
    if apply_noise:
        noise_mean=noise.mean(axis=0).mean(axis=0)
        noise_std=noise.std(axis=0).std(axis=0)
        noise_random=np.array([[np.random.normal(noise_mean[i][j].real,noise_std[i][j].real,msg_send_mod.shape[-1])+1j*np.random.normal(noise_mean[i][j].imag,noise_std[i][j].imag,msg_send_mod.shape[-1]) for i in range(noise.shape[-2])] for j in range(noise.shape[-1])])
    msg_recv_mod=np.array([bws[:,:,i].T@(msg_send_mod[i,:,:]+noise_random[i,:,:]) for i in range(bws.shape[2])]).astype('complex128')
    msg_recv_decode=np.array([decoder(bws,msg_recv_mod[:,:,i],apply_noise*noise.mean(axis=0).mean(axis=0)) for i in range(msg_recv_mod.shape[-1])])
    msg_recv=Modulation_type.demodulate(msg_recv_decode.reshape(-1), 'hard').reshape(msg_send.shape)
    # print(f'msg_send: {msg_send.shape}\nmsg_send_mod: {msg_send_mod.shape}\nnoise_random: {noise_random.shape}\nbws: {bws.shape}\nmsg_recv_mod: {msg_recv_mod.shape}\nmsg_recv_decode: {msg_recv_decode.shape}\nmsg_recv: {msg_recv.shape}')
    return msg_send,msg_recv


def cal_BER(msg_send,msg_recv):
    bit_num=msg_send.shape[0]*msg_send.shape[1]*msg_send.shape[2]
    bit_err=np.sum(msg_send!=msg_recv)
    BER=bit_err/bit_num
    print(f'Bit Errors: {bit_err}/{bit_num}={BER}')
    return BER

if __name__ == '__main__':
    np.random.seed(0)
    # generate_samples('./traces/ArgosCSI-8x6-2015-11-28-14-39-47_uhf_static/ArgosCSI-8x6-2015-11-28-14-39-47_static.hdf5',20)
    h5log, userCSI, noise, num_users, samps_per_user, timestep=openLog('./sample_traces/sample_20.hdf5',1000)
    bws=calc_martix(userCSI,0,10,'zf')
    msg_send,msg_recv=simulate(CDecoder, 64,bws,noise,1000,0.5,False)
    ber=cal_BER(msg_send,msg_recv)


