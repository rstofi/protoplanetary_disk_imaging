'''
 I, Kristóf Rozgonyi (rstofi@gmail.com) own the code under the The MIT License (MIT).
'''
import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.axes_grid1 import make_axes_locatable;

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
		
		im = plt.imshow(plot_data, cmap ='coolwarm', interpolation='nearest'); #coolwarm just a preference...

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
		
		cbar.set_label(r'$\Sigma$ [simulation units]', fontsize = 24);
		cbar.ax.tick_params(labelsize=18,pad=15);
		
		if output == None:
			### display image ###
			plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();
		
	def polar_plot(self,output=None,g=None,norm=None):
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
		
		### gamma correction ###
		if g != None:
			plot_data = np.power(plot_data,g);
		else:
			pass;
		
		### create grid ###
		r = np.arange(0,self.r)+self.r*0.1; #the center blank zone is 10% of the maximal r
		phi = np.linspace(0,2*np.pi,num=self.phi);
		
		phi_grid, r_grid = np.meshgrid(phi, r);
		x, y = r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid);
		
		### plot ###
		plt.figure(figsize=(10,10));
		im = plt.pcolormesh(x, y, plot_data, cmap ='coolwarm');
		
		ax = plt.gca();
		ax.tick_params(direction = 'inout', length= 6, width = 1.5, axis='both', which='major', labelsize=18, pad=15); #just a preference...
		
		plt.xlabel(r'X [simulation units]', fontsize = 24);
		plt.ylabel(r'Y [simulation units]', fontsize = 24);
		
		plt.xlim([-x.shape[0]*1.1,x.shape[0]*1.1]);
		plt.ylim([-x.shape[0]*1.1,x.shape[0]*1.1]);
		
		### set colorbar ###
		cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.1); #pad = 0 removes the gap between the plot and the colorbar
		
		cbar = plt.colorbar(im, cax=cax);
		
		cbar.set_label(r'$\Sigma$ [simulation units]', fontsize = 24);
		cbar.ax.tick_params(labelsize=18,pad=15);
		
		if output == None:
			### display image ###
			plt.tight_layout();
				
			plt.show();
			plt.close();
		else:
			plt.savefig(output, bbox_inches='tight');
			plt.close();
