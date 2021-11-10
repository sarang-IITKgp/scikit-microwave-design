This module contains customized functions for generating `Matplotlib` plots. 


###### polar_plot_colored_line(theta,R,T,ax, color_axis=None, 		use_colormap='jet',linewidth=3,title_name=None):
	Plots the line (theta,R) on polar ploton the axis 'ax', with continuous 
	color variation according to values in T and colormap as use_colormap


###### plot_colored_line(x,y,T,ax, color_axis=None, use_colormap='jet',linewidth=3,title_name=None,vmin=None,vmax=None):
    Plots the line (x,y), with continuous color variation on the axis 'ax' according to values in T and colormap as use_colormap.
    
###### plot_smith_chart(freq,S_p,fig,linewidth=5,use_colormap='hsv'):

	Generates a Smith chart and plots the values of S_p vs freq on the complex plane. 
	In general `freq` is supposed to be frequency, but it can be any other variable as well. 
	

	
