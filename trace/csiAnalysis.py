# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 09:45:06 2018

@author: jianding
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 09:09:37 2016

@author: Jian Ding
"""

import h5py
#from math import log, ceil, floor
import numpy as np
#import pandas as pd
import matplotlib
#matplotlib.use('WX') #why isn't this working? #matplotlib.rcsetup.all_backends #matplotlib.get_backend()
import matplotlib.pyplot as plt
#from warpNode import txgain2db
import lts
import time
import os

#import sys
#import os
#import glob

def samps2csi(samps, num_users,samps_per_user=224, fft_size=64):
	offset = 15+32 #todo: use findlts
	chunkstart = time.time()	
	#num_users = samps.shape[2]//samps_per_user
	print(samps.shape)
	
	# number of frames, number of antennas, number of samples in each frame, number of users, IQ
	usersamps = np.reshape(samps, (samps.shape[0],samps.shape[1],num_users,samps_per_user, 2) )	
	iq = np.empty((samps.shape[0],samps.shape[1],num_users,2,fft_size),dtype='complex64')

	for i in range(2): #2 seperate estimates
		iq[:,:,:,i,:] = (usersamps[:,:,:,offset+i*fft_size:offset+(i+1)*fft_size,0]+usersamps[:,:,:,offset+i*fft_size:offset+(i+1)*fft_size,1]*1j)*2**-15
	
	iq = iq.swapaxes(1,2)
	iq = iq.swapaxes(2,3)
	fftstart = time.time()
	
	csi = np.fft.fftshift(np.fft.fft(iq, fft_size, 4),4)*lts.lts_freq#*signal_max		
	endtime = time.time()
	print("chunk time: %f fft time: %f" % (fftstart - chunkstart, endtime -fftstart) )
	csi = np.delete(csi,[0,1,2,3,4,5,32,59,60,61,62,63],4) #remove zero subcarriers
	return csi	

	 

if __name__ == '__main__':
	starttime = time.time()	
	filename='/home/steven/Course/Network/dataset/traces/ArgosCSI-96x2-2016-03-31-15-53-27_Jian_left_to_right/ArgosCSI-96x2-2016-03-31-15-53-27_Jian_left_to_right.hdf5'
		
	if filename == '':
		filenames = [name for name in os.listdir('.') if name[-5:]=='.hdf5']
		filename = max(filenames)
	f = h5py.File(filename,'r')
	frame = 400 #1000 #660 #2599
	conjdata =[]
	zfdata =[]
	_here = os.path.dirname(__file__)
	csi = np.memmap(os.path.join(_here,'temp.mymemmap'), dtype='complex64', mode='w+', shape=(f['Pilot_Samples'].shape[0], f.attrs['num_mob_ant']+1, f['Pilot_Samples'].shape[1],52))
	csi[0:1000] = np.mean(samps2csi(f['Pilot_Samples'][:1000,:,:,:], f.attrs['num_mob_ant']+1, samps_per_user=f.attrs['samples_per_user']),2)
	
	chunk_num = int(csi.shape[0]/1000)
	
	print( "total:", chunk_num)
	for i in range(1, chunk_num):
		print("current:", i)
		csi[i*1000:i*1000+1000] = np.mean(samps2csi(f['Pilot_Samples'][i*1000:(i*1000+1000),:,:,:], f.attrs['num_mob_ant']+1, samps_per_user=f.attrs['samples_per_user']),2)
	#dimensions: frames, users, antennas, subcarriers
	csi[chunk_num*1000:] = np.mean(samps2csi(f['Pilot_Samples'][chunk_num*1000:csi.shape[0],:,:,:], f.attrs['num_mob_ant']+1, samps_per_user=f.attrs['samples_per_user']),2)

	
	timestep = f.attrs['frame_length']/20e6	
	f.close()
	# noise is saved as the last user
	noise = csi[:,-1,:,:] 
	userCSI = csi[:,:-1,:,:]
	

	i = 1
	plt.figure(1,figsize=(6,6))
	plt.clf()
	plt.plot(np.arange(0,userCSI.shape[0])*timestep,(np.angle(userCSI[:,0,i,25])))
	plt.grid('on')
	plt.xlabel('Time (s)')
	plt.ylabel('Angle')
	plt.title('Channel between base station antenna %d and user, subcarrier 1'%i)
	
	plt.figure(2,figsize=(6,6))
	plt.clf()
	plt.plot(np.arange(0,userCSI.shape[0])*timestep,np.absolute(userCSI[:,0,i,25]))
	plt.grid('on')
	plt.xlabel('Time (s)')
	plt.ylabel('Amplitude')
	plt.title('Channel between base station antenna %d and user, subcarrier 0'%i)
	plt.show()
	
	
	


	
 	