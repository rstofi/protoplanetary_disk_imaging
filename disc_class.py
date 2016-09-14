import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.axes_grid1 import make_axes_locatable;
import matplotlib.gridspec as gridspec;

def radial_pixel_to_R(r_arg,r_grid=256.0,r_min_val = 3.0,r_max_val = 50.0):
	"""
	Converts the simulation radila units to phisical dimension data.
	This method using logaritmic binning radial profile
	"""	
	r_min_val = r_min_val;
	r_max_val = r_max_val;
	
	rt = r_min_val*np.exp((r_arg/r_grid)*np.log(r_max_val/r_min_val));
	rt_plus1 = r_min_val*np.exp(((r_arg+1)/r_grid)*np.log(r_max_val/r_min_val));
	
	R = (2.0/3.0)*(np.power(rt_plus1,3.0)-np.power(rt,3.0))/(np.power(rt_plus1,2.0)-np.power(rt,2.0));
	
	return R;
	
class disc():
	"""
	This is a class for 2D FARGO protoplanetary discs
	"""

	def __init__(self,r,phi,fname):
		"""
		The r and phi grid define the disc structure.
		fname is the FARGO binary output file name
		"""
		self.r = r;
		self.phi = phi;
		self.fname = fname;
		
	def read(self):
		"""
		Return the r-phi matrix of the disc
		"""
		dat_file = open (self.fname,'rb');
		disc_dat = np.fromfile(dat_file,dtype = np.float64,);
		dat_file.close();

		return np.reshape(disc_dat,(self.r,self.phi));

	def norm_read(self,norm=None):
		"""
		Return the normalized disc or just read the disc
		"""		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			normed_data = np.divide(data,norm_data);
		
			return normed_data;
		
		else:
			return data;
			
	def radial_max_index(self,norm=None):
		"""
		Return the index of maximum element of the azimutaly averaged radial profile
		If the absolute position of the maximum is needed use: r_max = np.unravel_index(np.argmax(data),data.shape);
		azimutal_max_index can be defined the same way if needed
		"""
		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();

		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			data = np.divide(data,norm_data);

		phi_avg = np.mean(data, axis=1);
	
		r_max = np.argmax(phi_avg);
		
		return r_max;
				
	def radial_ring_average(self,dr=10,norm=None):
		"""
		Cut the radial values around the azimutaly averaged r max index in dr range, and then return the averaged azimutal profile
		Using of the normalization is required to get reasonable output!
		"""
		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();

		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			data = np.divide(data,norm_data);
		
		### cut around maximum ###
		r_max = pdisc.radial_max_index(norm=norm);
		
		dat_ring = data[r_max-dr:r_max+dr+1,:];
		
		### average ###
		phi_frame = np.mean(dat_ring, axis=0);
		
		
		return phi_frame;
		
	def azimutal_average(self,norm=None):
		"""
		Average the frame azimutaly and return the radial profile
		"""
		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();

		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			data = np.divide(data,norm_data);
			
		r_avg = np.average(data,axis=1);
		
		return r_avg;
	
	def r_array(self,phi,norm=None):
		"""
		Return the phi^th radial grid cell array of the disc
		Index starts from zero!
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			data = np.divide(data,norm_data);
		
		r_array = data[:,phi];
		
		return r_array;	

	def phi_array(self,r,norm=None):
		"""
		Return the r^th azimutal grid cell array of the disc
		Index starts from zero!
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		data = pdisc.read();

		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			data = np.divide(data,norm_data);
		
		phi_array = data[r,:];
		
		return phi_array;
			
### plot functions ###	

	def radial_plot(self,output=None,g=None,norm=None):
		"""
		Plot the disc phi versus radius (the plot sets are just optional)
		output: name of the output file for example: "example_disc.pdf"
		g: gamma correction raise the data to the power of g
		norm: other disc object that used to normal the data
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);

		### shift maximum to pi ###
		max_index = np.int(np.round(self.phi/2.0));
		max_indexes = np.unravel_index(np.argmax(plot_data),plot_data.shape);
		phi_max_index = max_indexes[1];

		roll = max_index-phi_max_index;
		if roll >= 0:
			plot_data = np.roll(plot_data,roll,axis=1);
		else:
			plot_data = np.roll(plot_data,self.phi+roll,axis=1);
		
		### gamma correction ###
		if g != None:
			plot_data = np.power(plot_data,g);
		else:
			pass;
		
		### set plot parameters ###
		x_len = len(plot_data[1]);
		y_len = len(plot_data[2]);
		
		plot_coords = [0,x_len,0,y_len];
		
		fig = plt.figure(figsize=(self.phi*0.025, self.r*0.025));
		plt.xlim([0,self.phi]);
		plt.ylim([0,self.r]);
		
		#im = plt.imshow(plot_data, cmap ='coolwarm', interpolation='nearest', rasterized=True); #coolwarm just a preference...
		im = plt.pcolormesh(plot_data, cmap ='coolwarm', rasterized=True); #coolwarm just a preference...

		plt.xlabel(r'$\varphi$', fontsize = 24);
		plt.ylabel(r'r [simulation units]', fontsize = 24);
		
		ax = plt.gca();
		ax.tick_params(direction = 'inout', length= 6, width = 1.5, axis='x', which='major', labelsize=20, pad=15); #just a preference...
		ax.tick_params(direction = 'inout', length= 6, width = 1.5, axis='y', which='major', labelsize=18, pad=15); #just a preference...
		
		ticks=[0,self.phi/4,self.phi/2,3*self.phi/4,self.phi];
		plt.xticks(ticks);
		
		xlabels = [r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'];
		ax.set_xticklabels(xlabels);
		
		### set colorbar ###
		cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.1); #pad = 0 removes the gap between the plot and the colorbar
		
		cbar = plt.colorbar(im, cax=cax);
		
		if norm == None:
			cbar.set_label(r'$\Sigma$ [simulation units]', fontsize = 24);
		else:
			cbar.set_label(r'$\Sigma$ \ $\Sigma_0$', fontsize = 24);
			
		cbar.ax.tick_params(labelsize=18,pad=15);
		
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();
		
	def polar_plot(self,output=None,g=None,norm=None,polar_plot_radius=None,log_transform=False):
		"""
		Plot the disc in polar coordinates
		output: name of the output file for example: "example_disc.pdf"
		g: gamma correction raise the data to the power of g
		norm: other disc object that used to normal the data
		"""
		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);

		### shift maximum to pi ###
		max_index = np.int(np.round(self.phi/2.0));
		max_indexes = np.unravel_index(np.argmax(plot_data),plot_data.shape);
		phi_max_index = max_indexes[1];

		roll = max_index-phi_max_index;
		if roll >= 0:
			plot_data = np.roll(plot_data,roll,axis=1);
		else:
			plot_data = np.roll(plot_data,self.phi+roll,axis=1);
		
		### gamma correction ###
		if g != None:
			plot_data = np.power(plot_data,g);
		else:
			pass;

		### create grid ###
		if log_transform == False:
			r = np.arange(0,self.r)+self.r*0.1; #the center blank zone is 10% of the maximal r
			phi = np.linspace(0,2*np.pi,num=self.phi);
		
			phi_grid, r_grid = np.meshgrid(phi, r);
			x, y = r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid);
		elif log_transform == True:
			### create grid ###
			R_polar = np.zeros(self.r);
			
			### set real radial values ###		
			for i in range(0,R_polar.size):
				R_polar[i] = radial_pixel_to_R(i);
			
			r_polar = R_polar+R_polar*0.1; #the center blank zone is 10% of the maximal r
			
			#r_polar = np.arange(0,self.r)+self.r*0.1; #the center blank zone is 10% of the maximal r
			phi_polar = np.linspace(0,2*np.pi,num=self.phi);
			
			phi_grid, r_grid = np.meshgrid(phi_polar, r_polar);
			x, y = r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid);
		else:
			raise Exception('Invalid argument value!');
		
		### plot ###
		plt.figure(figsize=(12,10));
		im = plt.pcolormesh(x, y, plot_data, cmap ='coolwarm');
		
		ax = plt.gca();
		ax.tick_params(direction = 'inout', length= 6, width = 1.5, axis='both', which='major', labelsize=18, pad=15); #just a preference...
		
		ax.axis('equal');
		
		if log_transform == True:
			plt.xlabel(r'X [transformed simulation units]', fontsize = 24);
			plt.ylabel(r'Y [transformed simulation units]', fontsize = 24);
		else:
			plt.xlabel(r'X [simulation units]', fontsize = 24);
			plt.ylabel(r'Y [simulation units]', fontsize = 24);
		
		### set plot area ###
		if polar_plot_radius == None and log_transform == False:
			plt.xlim([-x.shape[0]*1.1,x.shape[0]*1.1]);
			plt.ylim([-x.shape[0]*1.1,x.shape[0]*1.1]);
		elif polar_plot_radius == None and log_transform == True:
			plt.xlim([-np.max(R_polar)*1.1,np.max(R_polar)*1.1]);
			plt.ylim([-np.max(R_polar)*1.1,np.max(R_polar)*1.1]);
			pass;
		else:
			plt.xlim([-polar_plot_radius,polar_plot_radius]);
			plt.ylim([-polar_plot_radius,polar_plot_radius]);
		
		### set colorbar ###
		cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.1); #pad = 0 removes the gap between the plot and the colorbar
		
		cbar = plt.colorbar(im, cax=cax);
		
		if norm == None:
			cbar.set_label(r'$\Sigma$ [simulation units]', fontsize = 24);
		else:
			cbar.set_label(r'$\Sigma$ \ $\Sigma_0$', fontsize = 24);

		cbar.ax.tick_params(labelsize=18,pad=15);
		
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();

	def azimutaly_averaged_radial_profile(self,norm=None,output=None,log_transform=False):
		"""
		Plot: the azimutaly averaged radial profile.
		Using of the normalization is required to get reasonable output!
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);
		
		### azimutaly average ###
		r_avg = pdisc.azimutal_average(norm=norm);
				
		### plot ###
		fig = plt.figure(figsize=(10,5));
		
		if log_transform == False:
			r = np.arange(0,self.r,1);
		elif log_transform == True:
			r = np.zeros(self.r);
			
			### set real radial values ###		
			for i in range(0,r.size):
				r[i] = radial_pixel_to_R(i);
		else:
			raise Exception('Invalid argument value!');
		
		plt.plot(r,r_avg, color="#000080", linewidth=4);
		plt.tick_params(axis='both', which='major', labelsize=18);
		plt.xlim([np.amin(r),np.amax(r)]);
		plt.ylim([np.amin(r_avg)-(0.1*np.average(r_avg)),np.amax(r_avg)+(0.1*np.average(r_avg))]);
		if log_transform == False:
			plt.xlabel(r'r [simulation units]', fontsize=18);
		else:
			plt.xlabel(r'r [transformed simulation units]', fontsize=18);
		if norm == None:
			plt.ylabel(r'<$\Sigma$>$_{\varphi}$', fontsize=22);
		else:
			plt.ylabel(r'<$\Sigma $ \ $\Sigma_0$>$_{\varphi}$', fontsize=22);
		plt.grid(True);
		
		### save/show ###
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();
	
	def azimutal_profile_around_maxima(self,dr=10,norm=None,output=None):
		"""
		Plot the azimutal profile of the radial values around the azimutaly averaged r max index in dr range
		Using of the normalization is required to get reasonable output!
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);
		
		### shift maximum to pi ###
		max_index = np.int(np.round(self.phi/2.0));
		max_indexes = np.unravel_index(np.argmax(plot_data),plot_data.shape);
		phi_max_index = max_indexes[1];

		roll = max_index-phi_max_index;
		if roll >= 0:
			plot_data = np.roll(plot_data,roll,axis=1);
		else:
			plot_data = np.roll(plot_data,self.phi+roll,axis=1);

		### radial average ###
		phi_ring_avg = pdisc.radial_ring_average(dr=dr,norm=norm);
		
		### shift the maximum to the center ###
		max_index = np.int(np.round(self.phi/2.0));
		phi_ring_max_index = np.argmax(phi_ring_avg);

		roll = max_index-phi_ring_max_index;
		if roll >= 0:
			phi_ring_avg = np.roll(phi_ring_avg,roll);
		else:
			phi_ring_avg = np.roll(phi_ring_avg,self.phi+roll);

		### plot ###
		fig = plt.figure(figsize=(10,5));
			
		phi = np.arange(0,self.phi,1);

		plt.plot(phi,phi_ring_avg, color="#000080", linewidth=4);
		plt.tick_params(axis='x', which='major', labelsize=20);
		plt.tick_params(axis='y', which='major', labelsize=18);
		plt.xlim([0,self.phi]);
		plt.ylim([np.amin(phi_ring_avg)-(0.1*np.average(phi_ring_avg)),np.amax(phi_ring_avg)+(0.1*np.average(phi_ring_avg))]);
		plt.xlabel(r'$\varphi$ [rad]', fontsize=18);
		if norm == None:
			plt.ylabel(r'<$\Sigma$>$_{dr=%s}$' %dr, fontsize=22);
		else:
			plt.ylabel(r'<$\Sigma$ \ $\Sigma_0$>$_{dr=%s}$' %dr, fontsize=22);
		plt.grid(True);

		### set radialn ticks to x label ###
		ticks=[0,self.phi/4,self.phi/2,3*self.phi/4,self.phi];
		xlabels = [r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'];

		plt.xticks(ticks,xlabels);

		### save/show ###
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();

	def radial_and_azimutal_profile(self,dr=10,norm=None,output=None,log_transform=False):
		"""
		Plot: 
			-the azimutaly averaged radial profile 
			-the azimutal profile of the radial values around the azimutaly averaged r max index in dr range
			
		Using of the normalization is required to get reasonable output!
		"""
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);
		
		### shift maximum to pi ###
		max_index = np.int(np.round(self.phi/2.0));
		max_indexes = np.unravel_index(np.argmax(plot_data),plot_data.shape);
		phi_max_index = max_indexes[1];

		roll = max_index-phi_max_index;
		if roll >= 0:
			plot_data = np.roll(plot_data,roll,axis=1);
		else:
			plot_data = np.roll(plot_data,self.phi+roll,axis=1);
					
		### radial an dazimutal average arrays ###
		r_avg = pdisc.azimutal_average(norm=norm);
		phi_ring_avg = pdisc.radial_ring_average(dr=dr,norm=norm);
		
		### shift the maximum to the center ###
		max_index = np.int(np.round(self.phi/2.0));
		phi_ring_max_index = np.argmax(phi_ring_avg);

		roll = max_index-phi_ring_max_index;
		if roll >= 0:
			phi_ring_avg = np.roll(phi_ring_avg,roll);
		else:
			phi_ring_avg = np.roll(phi_ring_avg,self.phi+roll);
				
		### plot ###
		fig = plt.figure(figsize=(20.2,5));
		
		gs = gridspec.GridSpec(1,2);
				
		### spaceing between plots ###
		fig.subplots_adjust(hspace=0.2);
		
		### first -radial- plot ###
		ax1 = plt.subplot(gs[0, 0]);

		if log_transform == False:
			r = np.arange(0,self.r,1);
		elif log_transform == True:
			r = np.zeros(self.r);
			
			### set real radial values ###		
			for i in range(0,r.size):
				r[i] = radial_pixel_to_R(i);
		else:
			raise Exception('Invalid argument value!');
					
		ax1.plot(r,r_avg, color="#000080", linewidth=4);
		ax1.tick_params(axis='both', which='major', labelsize=18);
		ax1.set_xlim([np.amin(r),np.amax(r)]);
		ax1.set_ylim([np.amin(r_avg)-(0.1*np.average(r_avg)),np.amax(r_avg)+(0.1*np.average(r_avg))]);
		if log_transform == False:
			ax1.set_xlabel(r'r [simulation units]', fontsize=18);
		else:
			ax1.set_xlabel(r'r [transformed simulation units]', fontsize=18);
		if norm == None:
			ax1.set_ylabel(r'<$\Sigma$>$_{\varphi}$', fontsize=22);
		else:
			ax1.set_ylabel(r'<$\Sigma $ \ $\Sigma_0$>$_{\varphi}$', fontsize=22);
		ax1.grid(True);
		
		### second - azimutal - plot ###
		ax2 = plt.subplot(gs[0, 1]);
		
		phi = np.arange(0,self.phi,1);
		ax2.plot(phi,phi_ring_avg, color="#000080", linewidth=4);
		ax2.tick_params(axis='x', which='major', labelsize=20);
		ax2.tick_params(axis='y', which='major', labelsize=18);
		ax2.set_xlim([0,self.phi]);
		ax2.set_ylim([np.amin(phi_ring_avg)-(0.1*np.average(phi_ring_avg)),np.amax(phi_ring_avg)+(0.1*np.average(phi_ring_avg))]);
		ax2.set_xlabel(r'$\varphi$ [rad]', fontsize=18);
		if norm == None:
			ax2.set_ylabel(r'<$\Sigma$>$_{dr=%s}$' %dr, fontsize=22);
		else:
			ax2.set_ylabel(r'<$\Sigma$ \ $\Sigma_0$>$_{dr=%s}$' %dr, fontsize=22);
		ax2.grid(True);
				
		### write radialn x label ###
		ticks=[0,self.phi/4,self.phi/2,3*self.phi/4,self.phi];
		ax2.set_xticks(ticks);
		
		xlabels = [r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'];
		ax2.set_xticklabels(xlabels);
		
		### save/show ###
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();
			
	def frame_plot(self,dr=10,norm=None,output=None,g=None,plot_title=None,polar_plot_radius=None, log_transform=False):
		"""
		Plot: 
			-the azimutaly averaged radial profile 
			-the azimutal profile of the radial values around the azimutaly averaged r max index in dr range
			-the polar plot
			
		Using of the normalization is required to get reasonable output!
		"""
		
		### read data ###
		pdisc = disc(self.r,self.phi,self.fname);
		
		plot_data = pdisc.read();
		
		### normalization ###
		if norm != None:
			norm_data = norm.read();
			
			plot_data = np.divide(plot_data,norm_data);
		
		### gamma correction ###
		if g != None:
			plot_data = np.power(plot_data,g);
		else:
			pass;
		
		### shift maximum to pi ###
		max_index = np.int(np.round(self.phi/2.0));
		max_indexes = np.unravel_index(np.argmax(plot_data),plot_data.shape);
		phi_max_index = max_indexes[1];

		roll = max_index-phi_max_index;
		if roll >= 0:
			plot_data = np.roll(plot_data,roll,axis=1);
		else:
			plot_data = np.roll(plot_data,self.phi+roll,axis=1);
			
		### create grid ###
		if log_transform == False:
			r = np.arange(0,self.r)+self.r*0.1; #the center blank zone is 10% of the maximal r
			phi = np.linspace(0,2*np.pi,num=self.phi);
		
			phi_grid, r_grid = np.meshgrid(phi, r);
			x, y = r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid);
		elif log_transform == True:
			### create grid ###
			R_polar = np.zeros(self.r);
			
			### set real radial values ###		
			for i in range(0,R_polar.size):
				R_polar[i] = radial_pixel_to_R(i);
			
			r_polar = R_polar+R_polar*0.1; #the center blank zone is 10% of the maximal r
			
			#r_polar = np.arange(0,self.r)+self.r*0.1; #the center blank zone is 10% of the maximal r
			phi_polar = np.linspace(0,2*np.pi,num=self.phi);
			
			phi_grid, r_grid = np.meshgrid(phi_polar, r_polar);
			x, y = r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid);
		else:
			raise Exception('Invalid argument value!');
					
		### radial an dazimutal datasets ###
		r_avg = pdisc.azimutal_average(norm=norm);
		phi_ring_avg = pdisc.radial_ring_average(dr=dr,norm=norm);
		
		### shift the maximum to the center ###
		max_index = np.int(np.round(self.phi/2.0));
		phi_ring_max_index = np.argmax(phi_ring_avg);

		roll = max_index-phi_ring_max_index;
		if roll >= 0:
			phi_ring_avg = np.roll(phi_ring_avg,roll);
		else:
			phi_ring_avg = np.roll(phi_ring_avg,self.phi+roll);
		
		### plot ###
		fig = plt.figure(figsize=(10,15));
		
		gs = gridspec.GridSpec(3,2);
		
		if plot_title != None:
			plt.suptitle(plot_title, fontsize=24);
			plt.subplots_adjust(top=0.94);
		
		### spaceing between plots ###
		fig.subplots_adjust(hspace=0.3,wspace=0.3);
		
		### first -radial- plot ###
		ax1 = plt.subplot(gs[0, 0]);

		if log_transform == False:
			r = np.arange(0,self.r,1);
		elif log_transform == True:
			r = np.zeros(self.r);
			
			### set real radial values ###		
			for i in range(0,r.size):
				r[i] = radial_pixel_to_R(i);
		else:
			raise Exception('Invalid argument value!');
					
		ax1.plot(r,r_avg, color="#000080", linewidth=4);
		ax1.tick_params(axis='both', which='major', labelsize=18);
		ax1.set_xlim([np.amin(r),np.amax(r)]);
		ax1.set_ylim([np.amin(r_avg)-(0.1*np.average(r_avg)),np.amax(r_avg)+(0.1*np.average(r_avg))]);
		if log_transform == False:
			ax1.set_xlabel(r'r [transformed simulation units]', fontsize=18);
		else:
			ax1.set_xlabel(r'r [transformed simulation units]', fontsize=18);
		if norm == None:
			ax1.set_ylabel(r'<$\Sigma$>$_{\varphi}$', fontsize=22);
		else:
			ax1.set_ylabel(r'<$\Sigma $ \ $\Sigma_0$>$_{\varphi}$', fontsize=22);
		ax1.grid(True);
		
		### second - azimutal - plot ###
		ax2 = fig.add_subplot(212);
		phi = np.arange(0,self.phi,1);
		ax2.plot(phi,phi_ring_avg, color="#000080", linewidth=4);
		ax2.tick_params(axis='x', which='major', labelsize=20);
		ax2.tick_params(axis='y', which='major', labelsize=18);
		ax2.set_xlim([0,self.phi]);
		ax2.set_ylim([np.amin(phi_ring_avg)-(0.1*np.average(phi_ring_avg)),np.amax(phi_ring_avg)+(0.1*np.average(phi_ring_avg))]);
		ax2.set_xlabel(r'$\varphi$ [rad]', fontsize=18);
		if norm == None:
			ax2.set_ylabel(r'<$\Sigma$>$_{dr=%s}$' %dr, fontsize=22);
		else:
			ax2.set_ylabel(r'<$\Sigma$ \ $\Sigma_0$>$_{dr=%s}$' %dr, fontsize=22);
		ax2.grid(True);
		
		### second - azimutal - plot ###
		ax2 = plt.subplot(gs[0, 1]);
		
		phi = np.arange(0,self.phi,1);
		ax2.plot(phi,phi_ring_avg, color="#000080", linewidth=4);
		ax2.tick_params(axis='x', which='major', labelsize=20);
		ax2.tick_params(axis='y', which='major', labelsize=18);
		ax2.set_xlim([0,self.phi]);
		ax2.set_ylim([np.amin(phi_ring_avg)-(0.1*np.average(phi_ring_avg)),np.amax(phi_ring_avg)+(0.1*np.average(phi_ring_avg))]);
		ax2.set_xlabel(r'$\varphi$ [rad]', fontsize=18);
		ax2.set_ylabel(r'<$\Sigma$ \ $\Sigma_0$>$_{dr=%s}$' %dr, fontsize=22);
		ax2.grid(True);
				
		### write radialn x label ###
		ticks=[0,self.phi/4,self.phi/2,3*self.phi/4,self.phi];
		ax2.set_xticks(ticks);
		
		xlabels = [r'$0$',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'];
		ax2.set_xticklabels(xlabels);
		
		### third - polar - plot ###
		ax3 = plt.subplot(gs[1:, :]);
				
		im = ax3.pcolormesh(x, y, plot_data, cmap ='coolwarm', rasterized=True);
		
		ax = plt.gca();
		ax.tick_params(direction = 'inout', length= 6, width = 1.5, axis='both', which='major', labelsize=18, pad=15); #just a preference...
		
		ax.axis('equal');
		
		if log_transform == False:
			plt.xlabel(r'X [simulation units]', fontsize = 18);
			plt.ylabel(r'Y [simulation units]', fontsize = 18);
		else:
			plt.xlabel(r'X [transformed simulation units]', fontsize = 18);
			plt.ylabel(r'Y [transformed simulation units]', fontsize = 18);
		
		### set plot area ###
		if polar_plot_radius == None and log_transform == False:
			plt.xlim([-x.shape[0]*1.1,x.shape[0]*1.1]);
			plt.ylim([-x.shape[0]*1.1,x.shape[0]*1.1]);
		elif polar_plot_radius == None and log_transform == True:
			plt.xlim([-np.max(R_polar)*1.1,np.max(R_polar)*1.1]);
			plt.ylim([-np.max(R_polar)*1.1,np.max(R_polar)*1.1]);
			pass;
		else:
			plt.xlim([-polar_plot_radius,polar_plot_radius]);
			plt.ylim([-polar_plot_radius,polar_plot_radius]);		
		### set colorbar ###
		cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.1); #pad = 0 removes the gap between the plot and the colorbar
		
		cbar = plt.colorbar(im, cax=cax);
		
		if norm == None:
			cbar.set_label(r'$\Sigma$ [simulation units]', fontsize = 24);
		else:
			cbar.set_label(r'$\Sigma$ \ $\Sigma_0$', fontsize = 24);

		cbar.ax.tick_params(labelsize=18,pad=15);
		
		### save/show ###
		if output == None:
			### display image ###
			#plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();

