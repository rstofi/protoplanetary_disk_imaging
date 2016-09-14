import disc_class as dc;

normdisc = dc.disc(256,512,'initial_disc.dat');
pdisc = dc.disc(256,512,'disc_t=1.dat');

dr = 5; #the radial window radius around maximum value
g = None; #gamma correction
#pp_radius = 300; #the polar plot radius

plot_title='Plot title'; #the title of the plot

#pdisc.radial_plot(norm=normdisc,output='radial_example.pdf');
#pdisc.polar_plot(norm=normdisc,output='polar_example.pdf');
#pdisc.azimutaly_averaged_radial_profile(norm=normdisc,output='radial_profile_example.pdf');
#pdisc.azimutal_profile_around_maxima(norm=normdisc,dr=5,output='azimutal_profile_example.pdf');
#pdisc.radial_and_azimutal_profile(norm=normdisc,dr=5,output='radial_and_azimutal_example.pdf');
#pdisc.frame_plot(dr=dr,norm=normdisc,g=g,plot_title=plot_title,output='frame_plot_example.pdf');

#pdisc.radial_plot(norm=normdisc);
#pdisc.polar_plot(norm=normdisc,log_transform=True);
#pdisc.azimutaly_averaged_radial_profile(norm=normdisc,log_transform=True);
#pdisc.azimutal_profile_around_maxima(norm=normdisc,dr=5);
#pdisc.radial_and_azimutal_profile(norm=normdisc,dr=5,log_transform=True);
pdisc.frame_plot(dr=5,norm=normdisc,g=1.5,plot_title='Example frame plot', log_transform=True,polar_plot_radius=30);
