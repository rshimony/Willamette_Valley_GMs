Willamette Valley (WV) Ground Motions Project Workflow:

**	All input and output files referred here can be found in the ZENODO companion publication to the paper "Effects of Willamette Valley Sedimentary Structure on Ground Motions". 

1)	Creating rfiles:

We create 3 rfiles (binary filies holding the model structural and geophysical properties, used to run SW4 simulations):
1.	New WV basin model structure embedded in the USGS CVM (Stephenson et al., 2017). The WV basin filling is single material (SM) (600 m/s throughout).
Make this for both SRV and EugSpe basin depths.
Script: make_rfile_basin_model_sm.py
Inputs:
-	USGS CVM binary files (/wv_project/rfiles/cvm_bin_files/)
-	Ifile: text file describing the structure – every grid point is assigned a basin depth value (/wv_project/rfiles/ifiles/ifile_cartesian_4rfile_cal_eugspe.txt or /wv_project/ifiles/ifile_cartesian_4rfile_cal_srv).
2.	New WV basin model structure embedded in the USGS CVM (Stephenson et al., 2017). The WV basin filling is a 1D velocity gradient.
Make this for both SRV and EugSpe basin depths.
Script: make_rfile_basin_model_1d.py
Inputs:
-	USGS CVM binary files (/wv_project/rfiles/cvm_bin_files/)
-	Ifile: text file describing the structure – every grid point is assigned a basin depth value (/wv_project/rfiles/ifiles/ifile_cartesian_4rfile_cal_eugspe.txt or /wv_project/ifiles/ifile_cartesian_4rfile_cal_srv).
3.	Only USGS CVM.
Script: make_rfile_usgscvm.py
Inputs:
-	USGS CVM binary files (/wv_project/rfiles/cvm_bin_files/)

2)	Simulations data:
All synthetic .sac files (waveforms) are in: /wv_project/synthetic_data/sac_files/
All synthetic images (ground velocity magnitude) are in: /wv_project/synthetic_data/images/

3)	Observed data:
Download recordings from the three EQs: Salem. Scotts Mills and Springfield.
Downloading script: get_wf_recordings.py
This script download data from a single event at a time. For each EQ the inputs should be changed accordingly.
This script downloads three types of recordings:
1.	Raw waveforms data
2.	Data after response correction (unfiltered)
3.	Filtered data (low-pass 1Hz)
Saves metadata and plot figures of all waveforms in separate directories.
Data location: /wv_project/observed_data/

4)	Waveform and Spectra Comparisons:
Three output types:
1.	Model comparison - Makes three components waveform and spectra comparison figures between ALL DATA (4 new WV models, USGS CVM and observed data) for each station in a single event.
Script: make_wf&spec_model_comp.py
*Need to run this for every event separately.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + metadata/snr_stations.csv
-	Observed waveform data:
/wv_project/observed_data/ + *event* + waveforms/observed_data_1hz_snr/
-	USGS CVM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /steph/
-	EugSpe 1D synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /eugspe_1d/
-	EugSpe SM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + / eugspe_sm /
-	SRV 1D synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + / srv_1d /
-	SRV SM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + / srv_sm /
Output data location: /wv_project/wf_spec_comparisons/wf_spec_model_comp/
* Used in supplementary material

2.	Station comparison -  Makes three components waveform and spectra comparison figures for ONE SPECIFIC new WV model (1 (out of 4) new WV model, USGS CVM and observed data) for each station in a single event. 
Script: make_wf&spec_station_comp.py
* This script takes ONE EVENT and ONE MODEL at a time.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + metadata/snr_stations.csv
-	Observed waveform data:
/wv_project/observed_data/ + *event* + waveforms/observed_data_1hz_snr/
-	USGS CVM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /steph/
-	New WV model synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + *model*/
Output data location: /wv_project/wf_spec_comparisons/ wf_spec_station_comp/
* Used in supplementary material

3.	Station comparison - Makes ONE (East) component waveform and spectra comparison figures for ONE SPECIFIC new WV model (1 (out of 4) new WV model, USGS CVM and observed data) for each station in a single event.
Script: make_wf&spec_station_E_comp.py
* This script takes ONE EVENT and ONE MODEL at a time.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + metadata/snr_stations.csv
-	Observed waveform data:
/wv_project/observed_data/ + *event* + waveforms/observed_data_1hz_snr/
-	USGS CVM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /steph/
-	New WV model synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + *model*/
Output data location: /wv_project/wf_spec_comparisons/ wf_spec_station_comp/
* Used in WF and Spectra comparison figures in publication.

5)	FAS Analysis:
Calculates Fourier Amplitude Spectra (FAS) binned to specific frequencies (0.2, 0.3, 0.4, 0.5, 0.6, 0.7 Hz) from waveform data (observed, new WV model synthetics, USGS CVM synthetics). Saving the calculated binned FAS in a flatfile.
Script: fas_analysis_flatfile.py
*Need to run this for every event separately.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + metadata/snr_stations.csv
-	Valley SNR stations data:
wv_project/observed_data/ + *event* + metadata/valley_snr_st.csv
-	Observed waveform data:
/wv_project/observed_data/ + *event* + waveforms/observed_data_1hz_snr/
-	USGS CVM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /steph/
-	New WV model synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + *model*/
Output data location: /wv_project/ fas_analysis/



6)	IMs Analysis:
Calculates intensity measures (IMs - PGV, Arias Intensity, Cross Correlation) from waveform data (observed, new WV model synthetics, USGS CVM synthetics). Saving the calculated IMs in a flatfile.
Script: ims_analysis_flatfile.py
* This script takes ONE EVENT and ONE MODEL at a time.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + /metadata/snr_stations.csv
-	Valley SNR stations data:
wv_project/observed_data/ + *event* + /metadata/valley_snr_st.csv
-	Observed waveform data:
/wv_project/observed_data/ + *event* + waveforms/observed_data_1hz_snr/
-	USGS CVM synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + /steph/
-	New WV model synthetics: 
/wv_project/synthetic_data/sac_files/+ *event* + *model*/
Output data location: /wv_project/ im_analysis/

7)	Residuals Analysis:
Creates RESIDUALS and RATIOs values flatfiles. (FAS and PGV)
Residual values (observed - synthetics) and residuals ratios (USGS CVM / WV model).
Script: residuals_analysis_flatfile.py
	*Need to run this for every event separately.
	Inputs:
-	FAS data flatfiles:
/wv_project/fas_analysis/ + *event* + /*
-	IMs data flatfiles:
/wv_project/im_analysis/ + *event* + /*
	Outputs:
-	Flatfiles with all residual and ratio values for all models (USGS CVM, EugSpe 1D, EugSpe SM, SRV 1D, SRV SM). Separate files for FULL and VALLEY stations inventories.
-	A combined flatfile with all residuals and ratios from ALL MODELS.
-	A “MELTED” flatfile through IMs.
Output data location: /residual_ratio_analysis/ + *event* + /

8)	Digitized Wells and Velocity Gradient:
Using digitized well log .csv files to calculate a linear velocity gradient.
*The first part of this script is needed for constructing the rfiles (velocity gradient).
Script: make_gradient_from_digi_wells.py
Inputs:
-	Digitized well logs .csv files:
/wv_project/wells_data/
	Outputs:
-	All wells concatenated flatfile.
/wv_project/wells_data/
-	Figure of stacked digitized wells with mean data profile and data fit. Part of Fig 2 in paper.
/wv_project/figures/digi_wells/
-	Separate figures of all digitized wells.
/wv_project/figures/digi_wells/
-	Subplot figure of all digitized wells plotted together.
/wv_project/figures/digi_wells/

9)	Make wells depth and location maps:
Script: make_wells_map.py
Inputs:
-	Depth to units (EugSpe and SRV) in all wells flatfile.
/wv_project/wells_data/depth2units_wells.csv
	Outputs:
-	EugSpe wells map figure
-	SRV wells map figure
/wv_project/figures/wells_maps/

10)	Wave propagation maps:
Making ground velocity maps taken every 0.5s of SW4 simulation.
These images show the wave propagation pattern in the simulation.
Script: make_wave_propagation_maps.py
	Inputs:
-	Base map to plot SW4 image data (ground velocity data) on top (made in this script).
/wv_project/figures/wave_propagation/map4waveprop_final.png
-	SW4 ground velocity images.
/wv_project/synthetic_data/images/+ *event* + *model* + freq8_mags/*
Outputs:
-	Converted .csv files with SW4 images data. These used for easier plotting. Deleted after used (after map figures made) due to storage shortage.
/wv_project/synthetic_data/images/+ *event* + *model* + freq8_csv_data/*
-	Figures of ground velocity maps, depicting wave propagation.
/wv_project/figures/wave_propagation/+ *event* + *model* + freq8_mags_figs/






11)	Models Cross Sections:
SW4 cross section images of the x axis (E-W) from the 4 new WV models.
Script: plot_model_crosssection.py
*Takes ONLY ONE MODEL at a time.
Inputs:
-	SW4 E-W Vs cross section images. These images are separate from the rest of the SW4 images and should be used ONLY FOR CSs.
/wv_project/synthetic_data/sw4_ims_for_cs/ims4cs_+*event*+*model*+/*
	Outputs:
-	Cross section figures. These grouped later to form Figure 3.
/wv_project/figures/models_cs/

12)	Boxplots Figures:
Make 3 types of boxplots:
1.	Boxplots of residual values distribution in the WV between all models (4 new WV models and USGS CVM) for EACH IM in a SINGLE EVENT (only Valley stations). Grouped later with residuals maps.
2.	Boxplots of residual ratio values distribution between Valley and Mountain (outside of valley) stations, for EACH IM, EACH MODEL in EACH EVENT.
3.	STACKED boxplot figure, that compares between all models and all events. This is FIGURE 4 in the publication.
Script: make_boxplots.py
Inputs:
-	Models’ residual and ratio data flatfiles.
/wv_project/residual_ratio_analysis/+*event*+/+*model*+_res_ff.csv
-	Concatenated ALL MODELS residual and ratio data flatfiles.
/wv_project/residual_ratio_analysis/+*event*+/all_model_val_ff.csv
-	“Melted” residual and ratio data flatfiles.
/wv_project/residual_ratio_analysis/+*event*+/melt_val_ff.csv
	Outputs:
-	Type 1 (residual values distribution in the WV)
/wv_project/figures/boxplots/res_val_dist/
-	Type 2 (residual ratio values distribution, Valley and Mountain stations)
/wv_project/figures/boxplots/ratio_inout_val_dist/
-	Type 3 (STACKED boxplot figure). Figure 4.
/wv_project/figures/boxplots/








13)	Intensity Measures (IMs) Residual and Ratio Maps:
Plotting these on top of a regional map. The final maps will show the residual or ratio data points (for each recording station used) and the basin depth of the model.
Script: plot_im_maps.py
*Need to run this for every event separately.
Inputs:
-	Event catalog file.
/wv_project/observed_data/event_catalog_valevents.csv
-	Ifiles. For plotting basin depth in maps.
/wv_project/ifiles/
-	Full station inventory residuals and ratios flatfiles. Used for ratio maps.
/wv_project/residual_ratio_analysis/+*event*+/+*model*+_res_ff.csv
-	Valley station inventory residuals and ratios flatfiles. Used for residual maps.
/wv_project/residual_ratio_analysis/+*event*+/+*model*+_res_ff_val.csv
	Outputs:
-	Residuals maps for EACH IM, EACH MODEL and EACH EVENT.
/wv_project/figures/ims_maps/residuals_maps/
-	Ratios maps for EACH IM, EACH MODEL and EACH EVENT.
/wv_project/figures/ims_maps/ratios_maps/

14)	Station location maps:
Plotting stations’ locations on top of a regional map. Basin depth plotted as well for reference of station location in respect to the WV.
Script: plot_station_loc_map.py
*This script has three parts:
1. Plotting maps of a single station. This part needs to run separately for each event. For each event, the stations that recorded it will be plotted, with the epicenter of the according event.
2. Plotting a map with stations used in paper figures 9, 10, 11. Figure 8.
3. Checking which stations recorded all three events and saving these stations' data to a separate 'all_eq' flatfile. Plot maps of these stations' location with the three epicenters and basin depth plotted as well.
Inputs:
-	SNR stations data: 
/wv_project/observed_data/ + *event* + /metadata/snr_stations.csv
*The second part (paper stations) use the Salem station inventory.
-	Event catalog file.
/wv_project/observed_data/event_catalog_valevents.csv
-	Ifiles. For plotting basin depth in maps.
/wv_project/ifiles/
	Outputs:
-	Single station location maps
/wv_project/figures/st_loc_maps/+*event*+/
-	Paper stations map
/wv_project/figures/st_loc_maps/
-	All EQs stations data flatfile.
/wv_project/observed_data/all_eq_st.csv
-	All EQs station’s location maps.
/wv_project/figures/st_loc_maps/all_eqs/

15)	Regional data map:
Plotting a regional map with data used in this work.
Figure 1.
Script: plot_regional_data_map.py
Inputs:
-	Event catalog file.
/wv_project/observed_data/event_catalog_valevents.csv
-	SNR stations data: 
/wv_project/observed_data/ + *event* + /metadata/snr_stations.csv
-	All EQs stations data flatfile.
/wv_project/observed_data/all_eq_st.csv
-	Depth to units (EugSpe and SRV) in all wells flatfile.
/wv_project/wells_data/depth2units_wells.csv
-	Ifiles. For plotting basin depth.
/wv_project/ifiles/
-	Depth interpretation from gravity data (McPhee et al., 2014).
wv_project/grav_interp_depths/zzthick4.axyz

16)	Z1.0 and Z2.5:
Makes a Z2.5 and Z1.0 surface data file from data contained in an rfile (SW4 simulations data file).
Script: rfile2z_surface.py
Makes 4 Z_Surface files:
1)	From the USGS CVM rfile.
2)	From the WV model (SRV 1D) rfile.
3)	USGS CVM Valley.
4)	WV model (SRV 1D) Valley.
*Valley Files - z depth values only for grid point within the WV polygon, while the rest of the domain (outside the WV) is 0.
Each Z_Surface file contains both Z1.0 and Z2.5 values.
Inputs:
-	Rfiles:
/wv_project/rfiles/usgs_cvm/stephenson_rfile2.rfile
/wv_project/rfiles/wv_basin_model_1d/srv/rfile_srv_1d.rfile
	Outputs:
-	USGS CVM z_surface file:
/wv_project/z_surfaces/z_surface_cvm
-	WV model z_surface file:
/wv_project/z_surfaces/z_surface_wv_model
-	Valley USGS CVM z_surface file:
/wv_project/z_surfaces/ z_surface_valley_cvm.xyz
-	Valley WV model z_surface file:
/wv_project/z_surfaces/ z_surface_valley_wv_model.xyz
-	Z_Surface map figures:
/wv_project/z_surfaces/

17)	Basin Amplification Factor (BAF) and Shaking Ratio Maps – OpenQuake:
Create shaking ratio maps and calculate BAF using the OpenQuake engine with Parker et al., 2020 GMPE.
Script: baf_openquake.py
*Detailed explanation of the process in the script.
Inputs:
-	Vs30 grid (text file). Required for openquake shaking calculation. This grid is resampled to fit the Z_Surface grid size.
/wv_project/vs30/vs30_dom.xygt'
-	WV model Z Surface file.
/wv_project/z_surfaces/z_surface_valley_wv_model.xyz
	Outputs:
-	OpenQuake data input file. Created for M8 Rrup=100km. Full domain.
*Also an input for calculating shaking files.
/wv_project/shakemaps/openquake/oq_data_m8_rrup100.xyz
-	OpenQuake data input file, for ONLY WV grid points.
*Also an input for calculating shaking files.
/wv_project/shakemaps/openquake/oq_data_m8_rrup100_valley.xyz

-	OpenQuake shaking file. Full domain. Z2.5=nan
/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_nan.xyz
-	OpenQuake shaking file. WV only. Z2.5=values
/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_wv_only_z2pt5_vals.xyz
-	OpenQuake shaking file. Full domain. Z2.5=nan outside WV. Z2.5=values inside WV. 
/wv_project/shakemaps/openquake/oq_data_shaking_m8_rrup100_full_domain_z2pt5_vals.xyz

*Shaking ratio is calculated between:
oq_data_shaking_m8_rrup100_full_domain_z2pt5_vals.xyz
/
oq_data_shaking_m8_rrup100_full_domain_z2pt5_nan.xyz

-	Shakemap figures:
1.	shakemap with Z2.5 values
2.	shakemap with Z2.5 = nan
3.	shaking ratio map (values/nans OR basin/no basin)
wv_project/figures/baf_shakemaps/openquake/

18)	Basin Amplification Factor (BAF) and Shaking Ratio Maps – SW4 Synthetics:
SW4 synthetics data Basin Amplification Factor (BAF) calculation and plotting.
Script: baf_sw4_synts_data.py
This script reads in PGV images from SW4 simulations (hmax files), calculates PGV ratios between the selected new WV model (SRV 1D) and the USGS CVM, calculate BAF value from these and plots map of PGV (shaking) ratio with BAF value.
Inputs:
-	Event catalog file.
/wv_project/observed_data/event_catalog_valevents.csv
-	Base map image for the SW4 images to be plotted on top.
/wv_project/figures/wave_propagation/map4waveprop_final.png
-	PGV SW4 image files (hmax). These are images taken at the last cycle of the simulation, showing the peak ground velocity at each grid point throughout the domain. In this script, only map view (z=0) images are used.
/wv_project/synthetic_data/images/+ *event* +_+*model*+_maxdudt/image.cycle=4014.z=0.hmaxdudt.sw4img
	Outputs:
-	Converted images to .csv files for shaking ratio data, calculated from the hmax images.
/wv_project/figures/baf_shakemaps/sw4_synts/+ *event* +_baf.csv

19)	Final Sublot Figures:
Making final figure subplots combination to use in publication and supplementary material.
Script: final_fig_subplots.py
Each cell in this script can run separately and creates a set of (or a single) figures. 
*More detailed description on each subplot in the script
Describing each cell:
•	All EQs waveform and spectra comparison:
Inputs:
-	All EQs stations flatfile (stations that recorded all 3 EQs).
/wv_project/observed_data/all_eq_st.csv 
-	Three component waveform and spectra comp figures from the 3 EQs. Change the model every run of this cell.
/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/+ *event* +/+*model*+/*.png
-	All EQs station location maps.
/wv_project/figures/st_loc_maps/all_eqs/*

	Outputs:
-	All EQs waveform and spectra comparison figs.
/wv_project/figures/final_subplots/all_eq_st_comp/+*model*+/

•	Combining SINGLE IM residual maps from ALL EQs to one figure, along with a boxplot comparison of this IM.
Inputs:
-	Residuals maps
/wv_project/figures/ims_maps/residuals_maps/*+event+*/*+IM+*_val*.png
-	Boxplots comparisons
/wv_project/figures/boxplots/res_val_dist/*+event+*/*+IM+*_val.png
	Outputs:
-	IMs maps with boxplot comparison
/wv_project/figures/boxplots/res_val_dist/*+event+*/

•	All IMs Residuals Maps Comparison Figures.
Combining all IMs RESIDUALS maps of every EQ and MODEL pair.
Inputs:
-	 Residuals maps
/wv_project/figures/ims_maps/residuals_maps/*+event+*/
Outputs:
-	All IM residual maps comparison maps. For every EQ and model pair.
/wv_project/figures/final_subplots/all_im_maps/

•	Grouping selected IMs maps (f=0.3,0.5,0.7Hz,PGV) from the SRV 1D model simulations for all events. 
FIGURE 6 in publication.
Inputs:
-	Residual maps.
/wv_project/figures/ims_maps/residuals_maps/*+event+*/
	Output:
-	/wv_project/figures/final_subplots/grouped_res_maps/

•	All IMs Ratio Maps Comparison Figures.
Combining all IMs RATIOS maps of every EQ and MODEL pair.
Inputs:
-	 Ratio maps
/wv_project/figures/ims_maps/ratios_maps/*+event+*/
	Outputs:
-	All IM ratio maps comparison maps. For every EQ and model pair.
/wv_project/figures/final_subplots/all_im_maps/

•	Combining all 3-component waveform and spectra model comparison (from each station and each event) with its station location map.
Inputs:
-	Waveform and spectra comparison figures.
/wv_project/wf_spec_comparisons/wf_spec_model_comp/*+event+*/
-	Single station location maps.
/wv_project/figures/st_loc_maps/*+event+*/
	Outputs:
-	Waveform and spectra comparison with station location figures.
/wv_project/figures/final_subplots/wf_spec_comp_st_map/*+event+*/

•	Combining all ratio maps with Valley or Mountain stations boxplots.
The second part of this cell takes the (previously made, in this section) PGV, SRV 1D model combined boxplot ratio maps from the 3 EQs and grouping them together.
This is FIGURE 5 in the publication.
Inputs:
-	Residual ratio maps.
/wv_project/figures/ims_maps/ratios_maps/*+event+*/*ratio.png
-	Ratio in/out valley boxplots.
/wv_project/figures/boxplots/ratio_inout_val_dist/*+event+*/
	Outputs:
-	Ratio maps with Valley or Mountain stations boxplots.
/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/*+event+*/
	FIGURE 5:
		Inputs:
o	Salem PGV SRV 1D Ratio map + boxplot:
/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/salem/pgv_srv_1d_ratio.jpg
o	Scott Mills PGV SRV 1D Ratio map + boxplot:
/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/scottsmills/pgv_srv_1d_ratio.jpg
o	Springfield PGV SRV 1D Ratio map + boxplot:
/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/springfield/pgv_srv_1d_ratio.jpg
		Output:
o	/wv_project/figures/final_subplots/ratio_maps_inout_val_dist/

•	Cross Sections Figure. Making a grouped cross section figure, showing 3 selected W-E cross sections from all 4 new WV models. 
This is FIGURE 3 in the publication.
Inputs:
-	Model cross section figures.
/wv_project/figures/models_cs/*+model+*/
	Output:
-	/wv_project/figures/final_subplots/grouped_model_cs/

•	Wells data figure. Plotting 2 wells location maps (SRV and EugSpe penetrating wells) grouped with digitized well logs and averaged velocity gradient figure.
This is FIGURE 2 in the publication.
Inputs:
-	Wells maps.
/wv_project/figures/wells_maps/srv_wells_map.png
/wv_project/figures/wells_maps/eugspe_wells_map.png
-	Stacked digitized wells, with 1D averaged velocity gradient.
/wv_project/figures/digi_wells/digi_wells_stacked.jpg
	Outputs:
-	/wv_project/figures/final_subplots/wells_fig/wells_fig.png

•	Wave propagation subplots. Combining wave propagation maps from the selected model (SRV 1D) and the USGS CVM, from the 3 EQs. Making this combined figure for every imaged time step (0.5s) in the simulation. 
These figures are later made into wave propagation video (using FFMPEG) that is used in the supplementary materials.
Inputs:
-	Simulations ground velocity images (image for every 0.5s in the simulation).
/wv_project/figures/wave_propagation/*+event+*_*+model+*_freq8_mags_figs/*
	Outputs:
/wv_project/figures/final_subplots/wave_propagation/

•	Grouping two selected timesteps (14s, 25s) wave propagation maps (made in the section above).
This is FIGURE 7 in the publication.
Inputs:
-	Wave propagation subplots.
/wv_project/figures/final_subplots/wave_propagation/
	Outputs:
/wv_project/figures/final_subplots/wave_prop_14_25/wave_prop_14_25.jpg

•	Plotting three figures of selected waveform and spectra comparisons grouped.
This is FIGURES 9, 10, 11 in the publication.
Inputs:
-	Waveform and spectra comparison E component figures.
/wv_project/wf_spec_comparisons/wf_spec_station_E_comp/*+event+*/srv_1d/
	Outputs:
-	/wv_project/figures/final_subplots/paper_wf_spec_comp/

•	Making BAF summary figure. Grouping data maps used in the OpenQuake BAF process (Vs30 and Z2.5) with the OpenQuake BAF map (top row) with BAF maps from the simulations of the three EQs (Salme, Scotts Mills, Springfield) (bottom row).
This is FIGURE 12 in the publication.
Inputs:
-	Vs30 map.
/wv_project/vs30/vs30_map.png
-	Z2.5 map.
/wv_project/z_surfaces/z2pt5_valley_map.png
-	Shaking ratio – OpenQuake.
/wv_project/figures/baf_shakemaps/openquake/z2pt5_full_valley_vals_ratio.png
-	Shaking ratio maps – SW4.
/wv_project/figures/baf_shakemaps/sw4_synts/*+event+*_baf.png
	Output:
-	/wv_project/figures/final_subplots/BAF_fig/BAF_fig.png
![image](https://github.com/user-attachments/assets/7a0761a3-6d8e-4d8f-ac32-e0ada8b63fe5)
