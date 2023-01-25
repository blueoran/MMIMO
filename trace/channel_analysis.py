#******************************************************************************
# 
# File   : channel_analysis.py
# Authors:	Clay Shepard (cws [at] rice.edu)
# License:	Copyright 2015, Rice Efficient Computing Group. All rights reserved.
#
#******************************************************************************

import h5py
import numpy as np
import matplotlib.pyplot as plt
import lts
import time

"""Argos Channel Analysis Library"""

"""
	Notes:
		On large logs this takes a lot of memory -- try increasing swap (linux) or virtual memory (windows).
		In linux you can just create an empty file with dd and use mkswap.		
"""

def openLog(filename,chunk_size):
	"""Open an Argos HDF5 file and return important parameters for analysis. 
	

	Args: filename of Argos HDF5 log created using ArgosCC.

	Returns:
		h5log: Opened h5py file from log.
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier].
		noise: Complex numpy array with [Frame, BS Ant, Subcarrier]
		num_users: Number of users in trace
		samps_per_user: Number of raw IQ samples per user
		timestep: Time between CSI collection (frame length) in seconds
	"""
	h5log = h5py.File(filename,'r')
	samps_per_user = h5log.attrs['samples_per_user']
	num_users = h5log.attrs['num_mob_ant']
	timestep = h5log.attrs['frame_length']/20e6		
	csi,iq = samps2csi_large(h5log['Pilot_Samples'], num_users+1, samps_per_user,chunk_size=chunk_size)
	noise = csi[:,-1,:,:,:]
	userCSI = np.mean(csi[:,:num_users,:,:,:],2)
	return h5log, userCSI, noise, num_users, samps_per_user, timestep

def samps2csi_large(samps, num_users, samps_per_user=224, offset=47, chunk_size=1000):
	"""Wrapper function for samps2csi_main for to speed up large logs by leveraging data-locality. Chunk_size may need to be adjusted based on your computer."""
	
	if samps.shape[0]> chunk_size:
		#rather than memmap let's just increase swap... should be just as fast.
		#csi = np.memmap(os.path.join(_here,'temp1.mymemmap'), dtype='complex64', mode='w+', shape=(samps.shape[0], num_users, 2, samps.shape[1],52))
		#iq = np.memmap(os.path.join(_here,'temp2.mymemmap'), dtype='complex64', mode='w+', shape=(samps.shape[0], num_users, 2, samps.shape[1],64))
		csi = np.empty((samps.shape[0], num_users, 2, samps.shape[1],52), dtype='complex64')
		iq = np.empty((samps.shape[0], num_users, 2, samps.shape[1],64), dtype='complex64')
		chunk_num = samps.shape[0]//chunk_size
		for i in range(chunk_num):
			csi[i*chunk_size:i*chunk_size+chunk_size], iq[i*chunk_size:i*chunk_size+chunk_size] = samps2csi(samps[i*chunk_size:(i*chunk_size+chunk_size),:,:,:], num_users, samps_per_user=samps_per_user)
		csi[chunk_num*chunk_size:],iq[chunk_num*chunk_size:] = samps2csi(samps[chunk_num*chunk_size:,:,:,:], num_users, samps_per_user=samps_per_user)
	else:
		csi,iq = samps2csi(samps, num_users, samps_per_user=samps_per_user, offset=offset)
	return csi,iq
	
def samps2csi(samps, num_users, samps_per_user=224, offset=47):
	"""Convert an Argos HDF5 log file with raw IQ in to CSI. 
	
	Asumes 802.11 style LTS used for trace collection.
		
	Args:
		samps: The h5py or numpy array containing the raw IQ samples [Frame, BS Ant, Samples, I/Q].
		num_users: Number of users used in trace collection. (Last 'user' is noise.)
		samps_per_user: Number of samples allocated to each user in each frame.
		
	Returns:
		csi: Complex numpy array with [Frame, User, LTS 1/2, BS Ant, Subcarrier]. (Last 'user' is noise.)
		iq: Complex numpy array of raw IQ samples [Frame, User, LTS 1/2, BS Ant, samples]
	
	Example:
		h5log = h5py.File(filename,'r')
		csi,iq = samps2csi(h5log['Pilot_Samples'], h5log.attrs['num_mob_ant']+1, h5log.attrs['samples_per_user'])
	"""
	#could use findlts to find offset
	fft_size=64
	chunkstart = time.time()	
	#num_users = samps.shape[2]//samps_per_user
	usersamps = np.reshape(samps, (samps.shape[0],samps.shape[1],num_users,samps_per_user, 2) )
	
	iq = np.empty((samps.shape[0],samps.shape[1],num_users,2,fft_size),dtype='complex64')
	for i in range(2): #2 seperate estimates
		iq[:,:,:,i,:] = (usersamps[:,:,:,offset+i*fft_size:offset+(i+1)*fft_size,0]+usersamps[:,:,:,offset+i*fft_size:offset+(i+1)*fft_size,1]*1j)*2**-15
		
	iq = iq.swapaxes(1,2)  #this is trickery to keep numpy vectorized (fast), then keep the axis in the order we want
	iq = iq.swapaxes(2,3)
	fftstart = time.time()
	csi = np.fft.fftshift(np.fft.fft(iq, fft_size, 4),4)*lts.lts_freq
	endtime = time.time()
	print("Chunk time: %f fft time: %f" % (fftstart - chunkstart, endtime - fftstart) )
	csi = np.delete(csi,[0,1,2,3,4,5,32,59,60,61,62,63],4) #remove zero subcarriers
	return csi,iq	

def log2csi_hdf5(filename):
	"""Convert raw IQ log to CSI.
	
	Converts input Argos HDF5 trace to frequency domain CSI and writes it to the same filename with -csi appended.
	
	Args: filename of log
	
	Returns: None
	
	"""
	h5log = h5py.File(filename,'r')
	samps_per_user = h5log.attrs['samples_per_user']
	num_users = h5log.attrs['num_mob_ant']	
	
	#compute CSI for each user and get a nice numpy array
	#Returns csi with Frame, User, LTS (there are 2), BS ant, Subcarrier  
	#also, iq samples nic(Last 'user' is noise.)ely chunked out, same dims, but subcarrier is sample.
	csi,iq = samps2csi(h5log['Pilot_Samples'], num_users+1, samps_per_user) #+1 for noise
	
	# create hdf5 file to dump csi to
	h5f = h5py.File(filename[:-5]+'-csi.hdf5', 'w')	
	h5f.create_dataset('csi', data=csi)
	#todo: copy gains
	#todo check if file exists os.path.isfile(fname) if f in glob.glob(f):
	for k in h5log.attrs:
		h5f.attrs[k] = h5log.attrs[k]
	h5f.close()	
	h5log.close()	
	
def calCond(userCSI):
	"""Calculate the standard matrix condition number.

	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		
	Returns:
		condNumber_ave: The average condition number across all users and subcarriers.
		condNumber: Numpy array of condition number [Frame, Subcarrier].	
	"""
	condNumber = np.empty((userCSI.shape[0],userCSI.shape[3]),dtype='float32')
	for sc in range(userCSI.shape[3]):
			condNumber[:,sc] = np.linalg.cond(userCSI[:,:,:,sc])
	condNumber_ave = np.average(condNumber)
	return condNumber_ave,condNumber

def calDemmel(userCSI):
	"""Calculate the Demmel condition number.

	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		
	Returns:
		demmelNumber_ave: The average condition number across all users and subcarriers.
		demmelNumber: Numpy array of condition number [Frame, Subcarrier].	
	"""
	demmelNumber = np.empty((userCSI.shape[0],userCSI.shape[3]),dtype='float32')
	for sc in range(userCSI.shape[3]):
		
		# covariance matrix
		cov = np.matmul(userCSI[:,:,:,sc],np.transpose(userCSI[:,:,:,sc],[0,2,1]).conj())
		eigenvalues = np.abs(np.linalg.eigvals(cov))	
		demmelNumber[:,sc] = np.sum(eigenvalues,axis=1)/np.min(eigenvalues,axis=1)
	demmelNumber_ave = np.average(demmelNumber)
	return demmelNumber_ave,demmelNumber

def find_bad_nodes(userCSI, thresh=.1, user=0):
	"""
	Find bad nodes based on the correlation between a reference node and every other node.
	If all nodes are 'bad' then the reference node is assumed to be bad, and a new reference node is chosen.
	userCSI should be a channel trace taken in a stable environment, after AGC has settled.
	E.g., use a 2 s trace and only pass the second half: 	find_bad_nodes(userCSI[userCSI.shape[0]//2:])	
	Return a list of the indexes of bad nodes.
	
	Note that by only correlating one node with another, we maximize decorrelation when there is a sync issue.
	(If we use more antennas for reference, then the impact of the time sync is less.
	If we use fewer, then the bad antennas, which all of the same phase shift, dominate correlation.)
	An important trick is to normalize amplitudes of the CSI, so that the phase shift always has the same impact,
	regardless of channel gains (e.g. if the reference node has low or high signal strength relative to the bad node
	it would decrease the impact of a time sync issue).  
	This also allows a fixed threshold to work well, as a single sample shift of 4 antennas out of 8 relatively consistently
	causes a decorrelation of .25.  If half the samples (worst case) this is a .125 deviation from the mean.
	
	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		thresh: Lower threshold means more sensitivity to deviation in correlation.
		Higher is more false negatives, lower is more false positives.
		user: Index of user to use for correlation (default is first user).
		
	Returns:
		bad_nodes: List of the indices of nodes determined to be bad.
	"""
	num_nodes = userCSI.shape[2]//4
	if num_nodes == 1:
		return [True]
	corr = [None] * num_nodes
	node_good = [True] * num_nodes
	ref = 0
	ref_good = False
	userCSI /= np.abs(userCSI) #normalize amplitude so that the phase shift always has the same effect on correlation (important!)	
	

	while not ref_good:
		for n in range(num_nodes):
			if n != ref:
				sl = list(range(ref*4,ref*4+4)) + list(range(n*4,n*4+4))
				corr_vec = np.conj(userCSI[0,:,sl,:]) #why does the sl slice do a transpose?
				c = userCSI[:,:,sl,:] #:userCSI.shape[0]//2 #find_bad_nodes(userCSI[:userCSI.shape[0]//2])
				corr[n] = calCorr(c, corr_vec)[0][:,user]
				v = np.max(np.abs(corr[n] - np.mean(corr[n])))
				if  v > thresh:
					node_good[n] = False
		if np.sum(node_good) > 1: #this ref antenna is fine (or happened to have the same timing errors as another node...)
			ref_good = True
		else:
			#print('Warning! Node %d chosen for reference appears to be bad! Trying another.' % ref)
			ref +=1
			if ref == num_nodes:
					#print('No good nodes found!  If there are only 2 nodes, perhaps just 1 is bad.')
					return [False] * num_nodes
			node_good = [True] * num_nodes
			
	bad_nodes = [i for i, x in enumerate(node_good) if not x]
	return bad_nodes

def calCorr(userCSI, corr_vec):
	"""Calculate the per-frame correlation with a given correlation vector (for all users provided).
	
	To sub-sample in time/user/subcarrier slice the input userCSI before passing it.
	Note, to slice, e.g. for just one user, dimensionality has to be maintained, so slice like userCSI[:,[2],:,:]

	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		corr_vec: Vector to correlate with [BS Ant, User, Subcarrier]
		
	Returns:
		corr_total: Average correlation across subcarriers [Frame, User]
		sig_sc: Correlation on every [Frame, User, Subcarrier]
	
	Example:
		corr_total,sig_sc = calCorr(userCSI,np.transpose(np.conj(userCSI[frame,:,:,:]),(1,0,2)))
	"""
	sig_intf = np.empty((userCSI.shape[0],userCSI.shape[1],userCSI.shape[1],userCSI.shape[3]),dtype='float32')

	for sc in range(userCSI.shape[3]):
		sig_intf[:,:,:,sc] = np.abs(np.dot(userCSI[:,:,:,sc],corr_vec[:,:,sc])) / np.dot( np.abs(userCSI[:,:,:,sc]), np.abs(corr_vec[:,:,sc]) ) #TODO: can we get rid of the for loop?

	sig_sc = np.diagonal(sig_intf,axis1=1,axis2=2)
	sig_sc = np.swapaxes(sig_sc,1,2)
	corr_total = np.mean(sig_sc,axis=2)
	
	return corr_total, sig_sc
	

def calCoherence(userCSI, interval=5, timestep=.01):
	"""Calculate the expected correlation based coherence time for each user and sub-carrier.
	
	To sub-sample in time/user/subcarrier slice the input userCSI before passing it.
	Note, to slice, e.g. for just one user, dimensionality has to be maintained, so slice like userCSI[:,[2],:,:]

	Note that this is very computationally expensive, so it is recommended you subsample and use short interval.	
	
	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		interval: Length of maximum delay to calculate coherence for.
		timestep: Time between frames in the provided userCSI.
	
	Returns:
		corr: Expected correlation at every delay: [Delay,Correlation,User,Subcarrier]
		
	Example:
		subsample = 5
		corr = calCoherence(userCSI[::subsample],interval=2,timestep=timestep*subsample)
		
	"""
	num_samps = int(interval//timestep)
	corr = np.empty((userCSI.shape[0]-num_samps,num_samps,userCSI.shape[1],userCSI.shape[3]),dtype='float32')
	conjcsi = np.conj(userCSI[:,:,:,:]) #np.transpose(np.conj(userCSI[:,:,:,:]),(0,2,1,3))
	for frame in range(userCSI.shape[0]-num_samps):	 #leave enough of the array
		for user in range(userCSI.shape[1]):
			for sc in range(userCSI.shape[3]):
				corr[frame,:,user,sc] = np.abs(np.dot(userCSI[frame:frame+num_samps,user,:,sc],conjcsi[frame,user,:,sc]))/np.dot(np.abs(userCSI[frame:frame+num_samps,user,:,sc]),np.abs(conjcsi[frame,user,:,sc]))	
	return corr


def calCapacity(userCSI, noise, beamweights, downlink=False):
	"""Calculate the capacity of a trace with static beamweights.
	
	Apply a set of beamweights to a set of wideband user channels and calculate the shannon capacity of the resulting channel for every Frame.
	
	Note that if the beamweights are calculated with a frame from the trace, that frame will have unrealistic capacity since it will correlate noise as signal.
	
	Args:
		userCSI: Complex numpy array with [Frame, User, BS Ant, Subcarrier]
		noise: Complex numpy array with [Frame, BS Ant, Subcarrier]
		beamweights: Set of beamweights to apply to userCSI [BS Ant, User, Subcarrier]
		downlink: (Boolean) Compute downlink capacity if True, else Uplink
		
	Returns:
		cap_total: Total capacity across all users averaged over subarriers in bps/hz [Frame]
		cap_u: Capacity per user across averaged over subcarriers in bps/hz [Frame, User]
		cap_sc: Capacity per user and subcarrier in bps/hz [Frame, User, Subcarrier]
		SINR: Signtal to interference and noise ratio for each frame user and subcarrier [Frame, User, Subcarrier]
		cap_su_sc: Single user (no interference) capacity per subcarrier in bps/hz  [Frame, User, Subcarrier]
		cap_su_u: Single user (no interference) capacity averaged over subcarriers in bps/hz [Frame, User]
		SNR: Signtal to noise ratio for each frame user and subcarrier [Frame, User, Subcarrier]
	"""
	noise_bs_sc = np.mean(np.mean(np.abs(noise),0),0)  #average over time and the two ltss
	sig_intf = np.empty((userCSI.shape[0],userCSI.shape[1],userCSI.shape[1],userCSI.shape[3]),dtype='float32')
	noise_sc_u = np.empty((userCSI.shape[1],userCSI.shape[3]),dtype='float32')
	for sc in range(userCSI.shape[3]):
		sig_intf[:,:,:,sc] = np.square(np.abs(np.dot(userCSI[:,:,:,sc],beamweights[:,:,sc])))  #TODO: can we get rid of the for loop?
		noise_sc_u[:,sc] = np.dot(np.square(noise_bs_sc[:,sc]),np.square(np.abs(beamweights[:,:,sc]))) #noise is uncorrelated, and all we have is average power here (Evan wants to do it per frame, but I think that's a bad idea)
	
	#noise_sc_u *= 4 #fudge factor since our noise doesn't include a lot of noise sources

	sig_sc = np.diagonal(sig_intf,axis1=1,axis2=2)
	sig_sc = np.swapaxes(sig_sc,1,2)
	intf_sc = np.sum(sig_intf,axis=1+int(downlink)) - sig_sc  
	SINR = sig_sc/(noise_sc_u+intf_sc)
		
	cap_sc = np.log2(1+SINR)
	cap_u = np.mean(cap_sc,axis=2)
	cap_total = np.sum(cap_u,axis=1)
	
	SNR = sig_sc/noise_sc_u	 
	cap_su_sc = np.log2(1+SNR)
	cap_su_u = np.mean(cap_su_sc,axis=2)
	
	return cap_total,cap_u,cap_sc,SINR,cap_su_sc,cap_su_u,SNR


def calContCapacity(csi, conj=True, downlink=False, offset=1):
	"""Calculate the capacity of a trace with continuous beamforming.	

	For every frame in a trace, calculate beamweights (either conjugate or ZF), 
	apply them to a set of wideband user channels either from the same frame or some constant offset (delay), 
	then calculate the shannon capacity of the resulting channel.
	
	The main difference in uplink/downlink is the source of interference (and power allocation).
	In uplink the intended user's interference is a result of every other user's signal passed through that user's beamweights.
	In downlink the inteded user's interference is a result of every other user's signal passed through their beamweights (applied to the intended user's channel).
	
	Note that every user has a full 802.11 LTS, which is a repitition of the same symbol.
	This method uses the first half of the LTS to make beamweights, then applies them to the second half.
	Otherwise, noise is correlated, resulting in inaccurate results.
	
	Args:
		csi: Full complex numpy array with separate LTSs and noise [Frame, User, BS Ant, Subcarrier] (noise is last user)
		conj: (Boolean) If True use conjugate beamforming, else use zeroforcing beamforming.
		downlink: (Boolean) Compute downlink capacity if True, else Uplink
		offset: Number of frames to delay beamweight application.
		
	Returns:
		cap_total: Total capacity across all users averaged over subarriers in bps/hz [Frame]
		cap_u: Capacity per user across averaged over subcarriers in bps/hz [Frame, User]
		cap_sc: Capacity per user and subcarrier in bps/hz [Frame, User, Subcarrier]
		SINR: Signtal to interference and noise ratio for each frame user and subcarrier [Frame, User, Subcarrier]
		cap_su_sc: Single user (no interference) capacity per subcarrier in bps/hz  [Frame, User, Subcarrier]
		cap_su_u: Single user (no interference) capacity averaged over subcarriers in bps/hz [Frame, User]
		SNR: Signtal to noise ratio for each frame user and subcarrier [Frame, User, Subcarrier]
	"""
	csi_sw = np.transpose(csi,(0,4,1,3,2)) #hack to avoid for loop (matmul requires last two axes to be matrix) #frame, sc, user, bsant, lts
	noise = csi_sw[:,:,-1,:,:]  #noise is last set of data. #frame, sc, bsant, lts
	userCSI_sw = csi_sw[:,:,:-1,:,0] #don't include noise, use first LTS for CSI #frame, sc, user, bsant, lts
	
	
	noise_sc_bs = np.mean(np.mean(np.abs(noise),3),0)  #average over time and the two ltss
	
	if conj:	
		'''Calculate weights as conjugate.'''
		beamweights = np.transpose(np.conj(csi_sw[:,:,:-1,:,1]),(0,1,3,2))
	else:
		'''Calculate weights using zeroforcing.'''
		beamweights = np.empty((userCSI_sw.shape[0],userCSI_sw.shape[1],userCSI_sw.shape[3],userCSI_sw.shape[2]),dtype='complex64')	
		for frame in range(userCSI_sw.shape[0]):		
			for sc in range(userCSI_sw.shape[1]):
				beamweights[frame,sc,:,:] = np.linalg.pinv( csi_sw[frame,sc,:-1,:,1] ) #* np.linalg.norm(csi[frame,:4,0,:,sc]) #either this, or the noise power has to be scaled back accordingly
	if offset > 0:
		beamweights = np.roll(beamweights, offset, axis=0) #delay offset samples
	
	sig_intf = np.square(np.abs(np.matmul(userCSI_sw[offset:,:,:,:],beamweights[offset:,:,:,:])))

	noise_sc_u = np.transpose(np.sum(np.square(noise_sc_bs)*np.square(np.abs(np.transpose(beamweights,(0,3,1,2)) )),3),(0,2,1))	
	noise_sc_u = noise_sc_u[offset:]
	#noise_sc_u *= 4 #fudge factor since our noise doesn't include a lot of noise sources.  this should probably be justified/measured or removed

	sig_sc = np.diagonal(sig_intf,axis1=2,axis2=3)
	intf_sc = np.sum(sig_intf,axis=2+int(downlink)) - sig_sc  #lazy hack -- just sum then subtract the intended signal.
	SINR = sig_sc/(noise_sc_u+intf_sc)
		
	cap_sc = np.log2(1+SINR)
	cap_u = np.mean(cap_sc,axis=1)
	cap_total = np.sum(cap_u,axis=1)
	
	SNR = sig_sc/noise_sc_u
	cap_su_sc = np.log2(1+SNR)
	cap_su_u = np.mean(cap_su_sc,axis=1)
	
	return cap_total,cap_u,cap_sc,SINR,cap_su_sc,cap_su_u,SNR


def calExpectedCapacity(csi, user=0, max_delay=100, conj=True, downlink=False):
	"""Calculate the expected capacity for beamweights calculated with delayed stale CSI.
	
	
	Args:
		csi: Full complex numpy array with separate LTSs and noise [Frame, User, BS Ant, Subcarrier] (noise is last user)
		user: Index of user to compute for (note that other users still affect capacity due to their interference)
		max_delay: Maximum delay (in frames) to delay the beamweight computation.
		conj: (Boolean) If True use conjugate beamforming, else use zeroforcing beamforming.
		downlink: (Boolean) Compute downlink capacity if True, else Uplink
		
	Returns:
		cap: Average capacity across all frames for a given delay (in frames) in bps/hz [Delay]
	"""
	cap = []
	for d in range(max_delay):
		#print([d,time.time()])
		delayed = calContCapacity(csi, conj=conj, downlink=downlink, offset=d)
		cap.append(np.mean(delayed[1][:,user]))

	return cap






# ********************* Example Code *********************

if __name__ == '__main__':
	starttime = time.time()	
	show_plots = False	
	zoom = 0  #samples to zoom in around frame (to look at local behavior), 0 to disable
	pl = 0
	static = h5py.File('/home/steven/Course/Network/dataset/traces/ArgosCSI-8x6-2015-11-28-14-39-47_uhf_static/ArgosCSI-8x6-2015-11-28-14-39-47_static.hdf5','r')
	env = h5py.File('/home/steven/Course/Network/dataset/traces/ArgosCSI-8x6-2015-11-28-14-43-08_uhf_env/ArgosCSI-8x6-2015-11-28-14-43-08_env.hdf5','r')
	mobile = h5py.File('/home/steven/Course/Network/dataset/traces/ArgosCSI-8x6-2015-11-28-14-47-06_uhf_mob/ArgosCSI-8x6-2015-11-28-14-47-06_mob.hdf5','r')
	# env = h5py.File('logs/ArgosCSI-8x6-2015-11-28-14-43-08_env.hdf5','r')
	# mobile = h5py.File('logs/ArgosCSI-8x6-2015-11-28-14-47-06_mob.hdf5','r')

	
	frame = 10 #frame to compute beamweights from
	conjdata =[]
	zfdata =[]

	for h5log in [static, env, mobile]:
		#read parameters for this measurement data
		samps_per_user = h5log.attrs['samples_per_user']
		num_users = h5log.attrs['num_mob_ant']
		timestep = h5log.attrs['frame_length']/20e6
		noise_meas_en = h5log.attrs.get('measured_noise', 1)

		#compute CSI for each user and get a nice numpy array
		csi,iq = samps2csi(h5log['Pilot_Samples'], num_users+noise_meas_en, samps_per_user) #Returns csi with Frame, User, LTS (there are 2), BS ant, Subcarrier  #also, iq samples nicely chunked out, same dims, but subcarrier is sample.
		if zoom > 0:  #zoom in too look at behavior around peak (and reduce processing time)
			csi = csi[frame-zoom:frame+zoom,:,:,:,:]
			frame = zoom  #recenter the plots (otherwise it errors)		
		noise = csi[:,-1,:,:,:]	 #noise is last set of data.
		userCSI = np.mean(csi[:,:num_users,:,:,:],2) #don't include noise, average over both LTSs

		#example lts find:
		user = 0
		#so, this is pretty ugly, but we want all the samples (not just those chunked from samps2csi), so we not only convert ints to the complex floats, but also have to figure out where to chunk the user from.
		lts_iq = h5log['Pilot_Samples'][frame,0,user*samps_per_user:(user+1)*samps_per_user,0]*1.+h5log['Pilot_Samples'][frame,0,user*samps_per_user:(user+1)*samps_per_user,1]*1j
		lts_iq /= 2**15		
		# print(lts.findLTS(lts_iq))
		# offset = lts.findLTS(lts_iq)+32	 #Andrew wrote this, but I don't really like the way he did the convolve method...  works well enough for high SNRs.
		offset = 32	 #Andrew wrote this, but I don't really like the way he did the convolve method...  works well enough for high SNRs.
		print("LTS offset for user %d, frame %d: %d" % (user, frame, offset) )

		#compute beamweights based on the specified frame.
		conjbws = np.transpose(np.conj(userCSI[frame,:,:,:]),(1,0,2))		
		zfbws = np.empty((userCSI.shape[2],userCSI.shape[1],userCSI.shape[3]),dtype='complex64')	
		for sc in range(userCSI.shape[3]):
			zfbws[:,:,sc] = np.linalg.pinv(userCSI[frame,:,:,sc])
		# print(zfbws-conjbws)

		# import pdb;pdb.set_trace()
		downlink = True
		#calculate capacity based on these weights
		#these return total capacity, per-user capacity, per-user/per-subcarrier capacity, SINR, single-user capacity(no inter-user interference), and SNR
		conj = calCapacity(userCSI, noise, conjbws, downlink=downlink) #conjcap_total,conjcap_u,conjcap_sc,conjSINR,conjcap_su_sc,conjcap_su_u,conjSNR
		zf = calCapacity(userCSI, noise, zfbws, downlink=downlink) #zfcap_total,zfcap_u,zfcap_sc,zfSINR,zfcap_su_sc,zfcap_su_u,zfSNR 
		
		#plot stuff
		if show_plots:
			#Multiuser Conjugate
			plt.figure(1000*pl, figsize=(50,10))
			plt.plot(np.arange(0,csi.shape[0]*timestep,timestep)[:csi.shape[0]],conj[1])
			#plt.ylim([0,2])		
			plt.xlabel('Time (s)')
			plt.ylabel('Per User Capacity Conj (bps/Hz)')
			plt.show()
			
			#Multiuser Zeroforcing
			plt.figure(1000*pl+1, figsize=(50,10))
			plt.plot(np.arange(0,csi.shape[0]*timestep,timestep)[:csi.shape[0]],zf[1])
			#plt.ylim([0,2])			
			plt.xlabel('Time (s)')
			plt.ylabel('Per User Capacity ZF (bps/Hz)')		
			plt.show()
			
			#Single user (but show all users)
			plt.figure(1000*pl+2, figsize=(50,10))
			plt.plot(np.arange(0,csi.shape[0]*timestep,timestep)[:csi.shape[0]],conj[-2])
			#plt.ylim([0,2])		
			plt.xlabel('Time (s)')
			plt.ylabel('SUBF Capacity Conj (bps/Hz)')
			plt.show()
			pl += 1
			
		#save for exporting to matlab (prettier plots)
		conjdata.append(conj)	
		zfdata.append(zf)		
		
		del csi,iq #free the memory
		
	endtime = time.time()
	print("Total time: %f" % (endtime-starttime) )
	'''
	import scipy.io
	data = dict(timestep=timestep)
	data.update(dict(static_zf_cap_total=zfdata[0][0], static_zf_cap_u=zfdata[0][1],static_conj_cap_total=conjdata[0][0], static_conj_cap_u=conjdata[0][1], static_conj_cap_su_u=conjdata[0][-2]))
	data.update(dict(env_zf_cap_total=zfdata[1][0], env_zf_cap_u=zfdata[1][1], env_conj_cap_total=conjdata[1][0], env_conj_cap_u=conjdata[1][1], env_conj_cap_su_u=conjdata[1][-2]))
	data.update(dict(mobile_zf_cap_total=zfdata[2][0], mobile_zf_cap_u=zfdata[2][1],mobile_conj_cap_total=conjdata[2][0], mobile_conj_cap_u=conjdata[2][1], mobile_conj_cap_su_u=conjdata[2][-2]))	
	#data = dict(timestep=timestep, static_zf_cap_total=zfdata[0][0], static_zf_cap_u=zfdata[0][1],static_conj_cap_total=conjdata[0][0], static_conj_cap_u=conjdata[0][1], env_zf_cap_total=zfdata[1][0], env_zf_cap_u=zfdata[1][1],env_conj_cap_total=conjdata[1][0], env_conj_cap_u=conjdata[1][1], mobile_zf_cap_total=zfdata[2][0], mobile_zf_cap_u=zfdata[2][1],mobile_conj_cap_total=conjdata[2][0], mobile_conj_cap_u=conjdata[2][1], static_conj_cap_su_u=conjdata[0][-2], env_conj_cap_su_u=conjdata[1][-2], mobile_conj_cap_su_u=conjdata[2][-2])
	scipy.io.savemat('logs/capacity-frame_%d.mat' % frame, data)
	'''
	
'''
%example matlab script for loading the saved file
load capacity-frame_500
%timestep = 0.035
plot(0:timestep:timestep*(length(env_conj_cap_u)-1),env_conj_cap_u)
plot(0:timestep:timestep*(length(mobile_conj_cap_u)-1),mobile_conj_cap_u)
xlim([0,600])
ylim([0,5])
ylabel('Per User Capacity Conj (bps/Hz)')
xlabel('Time (s)')

figure(4)
plot(0:timestep:timestep*(length(mobile_conj_cap_u)-1),mobile_conj_cap_u(:,2))
xlim([0,120])
ylim([0,5])
xlabel('Time (s)')
ylabel('User Capacity Conjugate (bps/Hz)')
print -clipboard -dmeta %windows only
'''

#import os
#import glob
'''  
	#Example for simply converting raw IQ to CSI.
	import glob
	logdir = "logs/uhf_wb_traces_vito/"
	filenames = glob.glob(logdir+"*.hdf5")
	#filenames = ('ChannelTracesVitosLand/ArgosCSI-8x5-2015-12-19-00-00-29_good_uhf_mobile_2directionalpolarized_1staticmobile_2mobile',
	#			'ChannelTracesVitosLand/ArgosCSI-8x4-2015-12-18-22-34-02_good_static_uhf_vito_alldirectional',
	#			'ChannelTracesVitosLand/ArgosCSI-8x4-2015-12-18-22-53-16_good_uhf_envmobility_vito.hdf5',)

	for filename in filenames:
		print(filename)
		log2csi_hdf5(filename)
'''
