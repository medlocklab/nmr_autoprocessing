b = load('\Users\Alex\Documents\NMR RSPA\October_Data\october_data_raw_data.mat')
spectra = csm_nmr_spectra(b.X,b.ppm)
a = csm_calibrate_nmr(spectra, 'TSP') %tsp adjusting
c = csm_normalise(a.output.calibrated_spectra,'peak') %normalizing to peak
model = csm_rspa(c.output.normalised_spectra)
figure = csm_gui_alignment(c.output.normalised_spectra, 'rspa_model', model)