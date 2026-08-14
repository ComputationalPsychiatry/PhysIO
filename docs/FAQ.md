Frequently Asked Questions (FAQ)
================================


## 1. What is the PhysIO Toolbox?

PhysIO is a toolbox for model-based physiological noise correction of fMRI data.

PhysIO stands for Physiological Input/Output toolbox, which summarizes its core purpose. A quote from our [paper](http://dx.doi.org/10.1016/j.jneumeth.2016.10.019>):


> In short, the toolbox transforms physiological input, i.e. peripheral recordings, into physiological output, i.e. regressors encoding components of physiological noise [...] A modular Matlab implementation supports command-line operation and is compatible with all major fMRI analysis packages via the export of regressor text-files. For the Statistical Parametric Mapping [SPM](<http://www.fil.ion.ucl.ac.uk/spm>) software package in particular, PhysIO features a full integration as a Batch Editor Tool, which allows user-friendly, GUI-based setup and inclusion into existing preprocessing and modeling pipelines.


## 2. How does PhysIO differ from other toolboxes for physiological noise correction for fMRI using peripheral recordings?

Citing from the introduction of our [paper](http://dx.doi.org/10.1016/j.jneumeth.2016.10.019>) again 

>
> ### Highlights ###
* A Toolbox to integrate preprocessing of physiological data and fMRI noise modeling.
* Robust preprocessing via iterative peak detection, shown for noisy data and patients.
* Flexible support of peripheral data formats and noise models (RETROICOR, RVHRCOR).
* Fully automated noise correction and performance assessment for group studies.
* Integration in fMRI pre-processing pipelines as SPM Toolbox (Batch Editor GUI).
>


## 3. How do I cite PhysIO?

The **core references for PhysIO** are: 

1. Kasper, L., Bollmann, S., Diaconescu, A.O., Hutton, C., Heinzle, J., Iglesias, 
S., Hauser, T.U., Sebold, M., Manjaly, Z.-M., Pruessmann, K.P., Stephan, K.E., 2017. 
*The PhysIO Toolbox for Modeling Physiological Noise in fMRI Data*. 
Journal of Neuroscience Methods 276, 56-72. https://doi.org/10.1016/j.jneumeth.2016.10.019
    - *main PhysIO Toolbox reference, also a good starting point to learn about more about the methods in PhysIO (see next question)*
2. Frässle, S., Aponte, E.A., Bollmann, S., Brodersen, K.H., Do, C.T., Harrison, O.K., Harrison, S.J., Heinzle, J., Iglesias, S., Kasper, L., Lomakina, E.I., Mathys, C., Müller-Schrader, M., Pereira, I., Petzschner, F.H., Raman, S., Schöbi, D., Toussaint, B., Weber, L.A., Yao, Y., Stephan, K.E., 2021. *TAPAS: an open-source software package for Translational Neuromodeling and Computational Psychiatry*. Frontiers in Psychiatry 12, 857. https://doi.org/10.3389/fpsyt.2021.680811
    - *main TAPAS software collection reference, see [main TAPAS README](https://github.com/translationalneuromodeling/tapas#readme) for more details on TAPAS itself*

Please cite these papers if you use PhysIO in your work. Here is a minimum example snippet: 

> The analysis was performed using the Matlab PhysIO Toolbox ([1], version x.y.z,
> open-source code available as part of the TAPAS software collection: [2], 
> <https://www.translationalneuromodeling.org/tapas>)

If you use respiratory volume per time (RVT) regressors or preprocess respiratory traces 
for RETROICOR, please also cite:

3. Harrison, S.J., Bianchi, S., Heinzle, J., Stephan, K.E., Iglesias, S., Kasper L., 2021.
A Hilbert-based method for processing respiratory timeseries.
NeuroImage, 117787. https://doi.org/10.1016/j.neuroimage.2021.117787
    - *superior RVT computation, preprocessing of respiratory traces*

A **standard comprehensive snippet to include** in your method section could look like the following, assuming you use our specific implementation of RETROICOR, which uses Fourier expansions of different order for the estimated phases of cardiac pulsation (3rd order), respiration (4th order) and cardio-­respiratory interactions (1st order) following (Harvey et al., 2008), and include respiratory volume per time (RVT) as well as heart-rate variability (HRV) regressors.

> Physiological noise correction was performed using the Matlab PhysIO Toolbox ([1], version x.y.z,
> open-source code available as part of the TAPAS software collection: [2], 
> <https://www.translationalneuromodeling.org/tapas>). A RETROICOR model [4,5]) was employed, using Fourier expansions of different order for the estimated phases of cardiac pulsation (3rd order),
> respiration (4th order) and cardio-­‐respiratory interactions (1st order) [6]. Furthermore, respiratory volume per time (RVT, [3,7]) and heart rate variability (HRV, [8]) were modeled. 

4. Glover, G.H., Li, T.Q. & Ress, D. Image-­‐based method for retrospective correction
of PhysIOlogical motion effects in fMRI: RETROICOR. Magn Reson Med 44, 162-­‐
7 (2000).

5. Hutton, C. et al. The impact of PhysIOlogical noise correction on fMRI at 7 T.
NeuroImage 57, 101-­‐112 (2011).

6. Harvey, A.K. et al. Brainstem functional magnetic resonance imaging:
Disentangling signal from PhysIOlogical noise. Journal of Magnetic Resonance
Imaging 28, 1337-­‐1344 (2008).

7. Birn, R.M., Smith, M.A., Jones, T.B., Bandettini, P.A., 2008. The respiration response
function: The temporal dynamics of fMRI s ignal fluctuations related to changes in
respiration. NeuroImage 40, 644–654. doi:10.1016/j.neuroimage.2007.11.059
PhysIO Toolbox | Citing this work 20

8. Chang, C., Cunningham, J.P., Glover, G.H., 2009. Influence of heart rate on the
BOLD signal: The cardiac response function. NeuroImage 44, 857–869.
doi:10.1016/j.neuroimage.2008.09.029


If you use **noise ROIs (aCompCor) or 12/24 regressor motion modeling**, also include the
respective references:

9. Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise
correction method (CompCor) for BOLD and perfusion based fMRI. NeuroImage 37,
90–101. doi:10.1016/j.neuroimage.2007.04.042
    - *aCompCor*

10. Siegel, J.S., Power, J.D., Dubis, J.W., Vogel, A.C., Church, J.A., Schlaggar, B.L.,
Petersen, S.E., 2014. Statistical improvements in functional magnetic resonance
imaging analyses produced by censoring high-motion data points. Hum. Brain Mapp.
35, 1981–1996. doi:10.1002/hbm.22307
    - *Motion Regressors*

## 4. Where do I find more documentation for PhysIO?

* The [paper](http://dx.doi.org/10.1016/j.jneumeth.2016.10.019) describing its structure, objective and modules
* [README.md](https://gitlab.ethz.ch/physio/physio-doc/blob/master/README.md) in the main folder when downloading
    - For help on installation and getting started
* Quickstart
    - PDF (or markdown .md file)
    - Tutorial matlab-scripts
* Reference Manual (for developers)


## 5. I am using FSL, AFNI, BrainVoyager, etc., for my fMRI analyses. Do I need SPM for PhysIO to work?

No, the basic functionality of PhysIO, i.e. creating nuisance regressors for your GLM analysis, is available in plain Matlab. The following extra functionality related to automatizing and assessing noise correction, require the installation of SPM:

- GUI (SPM Batch Editor)
- Pipeline dependencies (automatic input of realignment parameters, feed-in of multiple regressors file to GLM)
- Model assessment via F-tests and automatic F-map/tSNR report
- Noise-ROIs model (read-in of nifti files via SPM)


## 6. I am using device X for physiological recordings. Does PhysIO support the physiological logfile format Y?

Currently, PhysIO natively supports the following physiological logfile types:

- Brain Imaging Data Structure (BIDS)
    - [Standard for peripheral recordings](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html)
    - both raw physiological traces and pre-computed pulse events are
    supported
- BioPac formats
    - Biopac `.mat`-export
        - assuming the following variables (as columns): `data`, `isi`, `isi_units`, `labels`, `start_sample`, `units`
        - See `tapas_physio_read_physlogfiles_biopac_mat.m` for details
    - Biopac `.txt`-export
        - assuming the following 4 columns, with one sample per row: respiratory, skin conductance (GSR),  cardiac (PPG), and trigger signal (on/off)
- General Electric
- Philips SCANPHYSLOG files (`SCANPHYSLOG<DateTime>.log`; all versions from release 2.6 to 5.3)
- Siemens formats
    - Siemens VB (files `.ecg`, `.resp`, `.puls`)
    - Siemens VD/VE (files `*_ECG.log`, `*_RESP.log`, `*_PULS.log`)
        - including CMRR-derived multiband-files
    - Siemens Human Connectome Project log files (preprocessed 3 column files `*_Physio_log.txt`)
    
See [Read-In of Logfiles](MANUAL_PART_READIN) for a detailed description of the expected file formats.

Furthermore, physiological recordings can be entered via a *custom* data format, i.e., providing one text file per device. The files should contain one amplitude value per line. The corresponding sampling interval(s) are provided as a separate parameter in the toolbox.

If your favourite logfile format is not supported, please contact the developers. We try everything to accomodate the read-in flexibility of the toolbox to your needs.


## 7. I am running the toolbox for a lot of subjects / on a remote server without graphics. Can I somehow reproduce the output figures relevant to assess the data quality?

Yes you can, using the toolbox function `tapas_physio_review`. This function takes the physio-structure as an input argument, which is per default saved as `physio.mat` in the specified output folder of your batch job.


## 8. How do I interpret the various output plots of the toolbox?

Have a look at our publication: _The PhysIO Toolbox for Modeling Physiological Noise in fMRI Data_ (http://dx.doi.org/10.1016/j.jneumeth.2016.10.019)

The figures there give a good overview of the toolbox output figures, in particular:
- Fig. S1 (supplementary): Philips Scan Timing Sync from `gradient_log` (explanation of `thresh.zero`, `thresh.sli`, `thresh.vol`, `thresh.vol_spacing`
- Fig. 3: Diagnostic Raw Time Series (cardiac cycle length curve, respiration histogram)
- Fig. 8C: Single Subject F-contrast results (cardiac regressors)
- Fig. 9: Group results/typical activation sites for F-contrasts of RETROICOR regressors (cardiac/resp/interaction)


## 9. I want to access subject's physiological measures, e.g. heart rate or respiratory volume (per time), before they enter the regressors. Where can I do that?

All intermediate data processing steps (e.g. filtering, cropping) of the peripheral data, including the computation of physiologically meaningful time courses, such as heart rate and respiratory volume, are saved in the substructure `ons_secs` ("onsets in seconds) of the physio-structure mentioned in question 7. This structure is typically saved in a file `physio.mat`.

`physio.ons_secs` then contains the different time courses, cropped to the acquisition window synchronized to your fMRI scan (the same values before synchronization/cropping, is found in `physio.ons_secs.raw`). Here are the most important ones:

- `ons_secs.t`              	 = [];  % time vector corresponding to c and r
- `ons_secs.c`              	 = [];  % raw cardiac waveform (ECG or PPU)
- `ons_secs.r`              	 = [];  % raw respiration amplitude time course
- `ons_secs.cpulse`         	 = [];  % onset times of cardiac pulse events (e.g. R-peaks)
- `ons_secs.fr`                  = [];  % filtered respiration amplitude time series
- `ons_secs.c_sample_phase`      = [];  % phase in heart-cycle when each slice of each volume was acquired
- `ons_secs.r_sample_phase`      = [];  % phase in respiratory cycle when each slice of each volume was acquired
-   `ons_secs.hr`                  = [];  % [nScans,1] estimated heart rate at each scan
-   `ons_secs.rvt`                 = [];  % [nScans,1] estimated respiratory volume per time at each scan
-  `ons_secs.c_outliers_high`     = [];  % onset of too long heart beats
-  `ons_secs.c_outliers_low`      = [];  % onsets of too short heart beats 
-  `ons_secs.r_hist`              = [];  % histogram of breathing amplitudes
    
For a detailed list of all properties and their documentation, read the source code of `tapas_physio_new.m`


## 10. What is the order of the regressor columns in the multiple regressors file?

This depends on the physiological models (and their order) specified in the `model`-submodule of physio (or in the batch editor). The general order is outlined in Fig. 7A of the  [Main PhysIO Toolbox Paper](http://dx.doi.org/10.1016/j.jneumeth.2016.10.019). The []-brackets indicate the number of regressors:

1. RETROICOR cardiac regressors [2 x nOrderCardiac]
2. RETROICOR respiratory regressors [2 x nOrderRespiratory]
3. RETROICOR cardXResp interaction regressors [4 x nOrderCardiacXRespiratory]
4. HRV [nDelaysHRV]
5. RVT [nDelaysRVT]
6. Noise ROIs (PCA signatures and mean of each region) [nNoiseROIs x (nComponents+1)]
7. Other (included other text file) [nColumnsOtherFile]
8. Motion [6 or 12 or 24, depending on motion model]

If any of the models was not specified, the number of regressors is reduced accordingly.


## 11. How do I know whether the physiological noise correction worked?

The best way to assess the quality of the correction is an F-test over the respective physiological noise model regressors in the design matrix. Luckily, if you use SPM, the toolbox can create these contrasts and corresponding output plots with overlays of your brain automatically via calling the following function in the Matlab command window:

```
  args = tapas_physio_report_contrasts(...
                  'fileReport', 'physio.ps', ...
                  'fileSpm', 'analysisFolder/SPM.mat', ...
                  'filePhysIO', 'analysisFolder/physio.mat', ...
                  'fileStructural', 'anatomyFolder/warpedAnatomy.nii')
```

Of course, you will have to adapt all paths to your `SPM.mat`, `physio.mat` and `anatomy.nii` files. There are more parameters to set (e.g. F-contrast thresholds), type `help tapas_physio_report_contrasts` for a list of options.

There should be whole-brain multiple-comparison corrected "activation" in physiological noise sites (similar to Fig. 8C or 9 in our [paper](https://doi.org/10.1016/j.jneumeth.2016.10.019).

If your F-contrast results differ or are absent, have a look at the *Diagnostic raw physiological time series*-plot and check whether it resembles Fig. 3 in the paper or whether there are any suspicious spikes in the heart cycle length.

Other than that, scan timing synchronisation is a major source of error, so always check the *Cutout actual scans* plot, whether the curves and scan events, TR etc. make sense.


## 12. Philips: I would like to use the gradient log for timing synchronization, but how do I set the thresholds?

Have a look at the following figure:

![](images/FigureS1_GradientLogThresholds.png)

This figure can be found as [figure S1](https://www.sciencedirect.com/science/article/pii/S016502701630259X#sec0125) in the supplementary material of our [paper](https://doi.org/10.1016/j.jneumeth.2016.10.019).

The following heuristics might help with the threshold settings in the `sync`
structure:

1. Note that these thresholds have to be set correctly only once for each
functional sequence, i.e., usually once per study. Even small changes to scan
geometry (e.g. slice tilt) between subjects shouldn't affect them significantly.
2. Setting the thresholds is an iterative procedure. You might start with the
defaults, probably running into an error or warning (`Warning: Invalid
  MinPeakHeight. There are no data points greater than MinPeakHeight.` or `Not
  enough volume/slice scan events found`). Then you inspect the figure output
  resembling the one above and adjust (usually lower) the thresholds in the
  order mentioned below.
3. There are three time courses in the upper of the three subplots shown in the
figure. These time courses show the traces of the three gradient directions
`x`,`y`,`z`. Choose the one as `sync.grad_direction` parameter that has the
highest peaks and most regular features reflecting slice and/or volume scan
events.
4. `sync.zero` has to be smaller than `sync.slice` and `sync.vol`. It should
be about 4/5 of the typical peak height in the gradient trace. Note that you can
set this thresholds (and all other) either in absolute values or relative to the
maximum peak height. Set a value below 1, if you prefer the latter.
5. `sync.slice` should be about 9/10 of the typical peak height of a slice scan
event.
6. `sync.vol`, if you set it, should be larger than `sync.slice`. It should be
9/10 of the peak height that stands out at the beginning of a volume, and is
followed by some dozens of smaller peaks (for the slices) typically. It might
be, however, that there is no such peak marking the start of a volume. If so,
you might try `sync.vol_spacing` or leave it empty and rely on the slice
thresholds exclusively
6. `sync.vol_spacing`, if set, should reflect the temporal spacing (in seconds),
between the end of the previous volume and the start of the next one. The figure
above gives some idea how to do that based on the bottom subplot that shows peak
onset differences. If once every few seconds (your `TR`) you find an exposed
peak, its height will give you the value for `sync.vol_spacing` (maybe reduce
it by about 5-10ms to allow for timing inaccuracies).


## 13. How do I know which logfile type ('vendor') I have to choose?

- Typically, you will know your scanner manufacturer or the supplier of your peripheral recording device. The currently supported vendors can always be found in the SPM Batch Editor, as dropdown options for the vendor parameter in any PhysIO batch, and are also listed as cases in `tapas_physio_read_physlogfiles.m`.
- For Siemens, since there are a couple of formats, it is often helpful to check the extensions of the files (or the file name structure in general) see question 7.
- Sometimes you will have to look in the log files themselves and compare them to the examples provided on the [Data Section](https://www.tnu.ethz.ch/en/software/tapas/data.html) of our homepage.


## 14. What does Parameter *XY* mean and what is its best setting?

Before you ask us directly, there are two simple ways to find out more about the parameters and options of the PhysIO toolbox:

- In SPM, you can use the Batch Editor as a Help GUI directly. If you open or create a TAPAS PhysIO Batch and click on any parameter, there will be useful information about its meaning and suitable values in the help window, located in the lower part of the Batch Editor.
- Within Matlab, type `edit tapas_physio_new`. This constructor function lists all parameters of the physio-structure with inline comments on their purpose and possible values.

## 15. I don't have Matlab, can I still use the PhysIO Toolbox?

Yes, our friends at [NeuroDesk](https://neurodesk.github.io/) provide a docker container that contains the latest standalone version of SPM with inbuilt PhysIO. NeuroDesk runs in a browser and on any major operating system, and no Matlab license is required. For more details on PhysIO in NeuroDesk, check out their [Tutorials](https://neurodesk.github.io/tutorials/) on functional imaging.

Internally, this is using a compiled version of SPM with PhysIO. Contact us if you want to compile SPM with PhysIO yourself - for example, to deploy your own modifications of PhysIO to collaborators without Matlab.

Once SPM+PhysIO is compiled and you have the right Matlab Runtime Environment installed on your computer (determined by the Matlab version used for compilation), you can run a PhysIO SPM Batch file (`.m` or `.mat`) by calling

```
./run_spm12.sh ../MATLAB_Runtime/v99 batch physio_batch.mat
``` 

If you omit the batch file (but keep the keyword `batch`), this will just open the Batch Editor to interactively create and run a PhysIO batch.

## 16. PhysIO fails already when running your provided example data, e.g., with "Out of Memory" errors. What happened?

All PhysIO examples are thoroughly tested before each release using a unit/integration testing framework, but of course we cannot catch every bug in advance. 

The most likely explanation, however, is that the setup of PhysIO on your machine was not successful. Here are a few things you can try:

- Did you run `tapas_physio_init()` and did it give any warnings?
- Do you have the subfolders of SPM in your Matlab path = e.g., you used `addpath(genpath('path/to/spm12'))`? 
    - If yes, try `rmpath(genpath('path/to/spm12'));addpath('path/to/spm12')` and run the example again.
    - More details on this issue can be found in our [GitHub Forum](https://github.com/translationalneuromodeling/tapas/issues/166).


## 17. Which models do I have to include in my physiological regressor matrix? And which number of regressors (model order / delays) per model?

The question about how many and which regressors you need for successful noise removal depends a lot on your experimental design and your research questions. Here are a couple of rules of thumb:

1. To have sufficient degrees of freedom to estimate the parameters of the general linear model (GLM), it is often recommended to have a least 10 times as many datapoints as regressors in your model. That means, for 26 regressors in the GLM, it is advisable to have 260 volumes or more.
2. 26 regressors is the standard output if you include 
    1. RETROICOR (3rd order cardiac = 6 regressors (1 cosine and one sine per order), 4th order respiratory = 8 regressors, and 1st order multiplicative terms (4: cos*cos, sin*sin, cos*sin, sin*cos), 
    2. HRV (1 regressor), 
    3. RVT (1) and 
    4. motion realignment parameters (6)
    5. This list also gives an indication of the order of the regressors in the matrix, but see [question 10](https://tnurepository.ethz.ch/physio/physio-public/-/wikis/FAQ/edit#10-what-is-the-order-of-the-regressor-columns-in-the-multiple-regressors-file) for details.
3. You can now reduce the number of regressors by the following observations:
    1. The RETROICOR 3/4/1 order is taken from a paper that optimized physiological noise removal in the brainstem [1]. You might not need the full model, if you are interested in other brain areas
    2. For example, in our PhysIO paper, figure 9 [2], where we evaluated this model for 35 subjects, you   can see that the multiplicative terms explain noise mostly in the midbrain, acquaeduct and more inferior parts of the brainstem.
    3. Furthermore, also in that figure, you can see that the cardiac RETROICOR terms explain most of the variance (temporal SNR (tSNR) gains of up to 70%, compared to only motion correction), whereas the effect of both respiratory and multiplicative terms is one order of magnitude smaller (5% and 3% tSNR gains, respectively).
    4. So, you could probably leave out the multiplicative terms and reduce the number of regressors by 4, or reduced the respiratory terms to 2nd or 1st order, reducing them by 4 or 6, respectively. I have also seen 2nd cardiac order models, reducing by a further 2 regressors.
    5. In total, you would end up with 14 regressors for a RETROICOR 3/2/1, 10 regressors for RETROICOR 3/2/0 and 8 regressors for 2/2/0. In each case, you would still add the 8 regressors from HRV, RVT and motion.
4. To get a first idea which sets of regressors contribute to noise removal for your data, you can run F-Tests over the columns for the respective regressors only (see also [question 11](https://tnurepository.ethz.ch/physio/physio-public/-/wikis/FAQ/edit#11-how-do-i-know-whether-the-physiological-noise-correction-worked) and see where in the brain they explain significant variance (you can also use tapas_physio_compute_tsnr_gains to compute tSNR as in our figure).
    5. To formally compare whether it’s warranted to include a set of regressors into your GLM, you would have to do model comparison (as in [1]). There is a toolbox for SPM that can do this automatically for the kind of regressor matrices that are output by PhysIO [3].


[1] Harvey, A.K., Pattinson, K.T.S., Brooks, J.C.W., Mayhew, S.D., Jenkinson, M., Wise, R.G., 2008. Brainstem functional magnetic resonance imaging: Disentangling signal from physiological noise. Journal of Magnetic Resonance Imaging 28, 1337–1344. https://doi.org/10.1002/jmri.21623

[2] https://www.sciencedirect.com/science/article/pii/S016502701630259X#fig0045 from: Kasper, L., Bollmann, S., Diaconescu, A.O., Hutton, C., Heinzle, J., Iglesias, S., Hauser, T.U., Sebold, M., Manjaly, Z.-M., Pruessmann, K.P., Stephan, K.E., 2017. The PhysIO Toolbox for Modeling Physiological Noise in fMRI Data. Journal of Neuroscience Methods 276, 56–72. https://doi.org/10.1016/j.jneumeth.2016.10.019

[3] Soch, J., Allefeld, C., 2018. MACS – a new SPM toolbox for model assessment, comparison and selection. Journal of Neuroscience Methods 306, 19–31. https://doi.org/10.1016/j.jneumeth.2018.05.017


## 18. I cannot find the answer to my question in the FAQ. Whom do I ask for help?

We are very happy to provide support on how to use the PhysIO Toolbox. However, 
as every researcher, we only have a limited amount of time. So please excuse, if 
we might not provide a detailed answer to your request, but just some general 
pointers and templates. Before you contact us, please try the following:

1. A first look at the [FAQ](https://gitlab.ethz.ch/physio/physio-doc/wikis/FAQ) 
   (which is frequently extended) might already answer your questions.
2. For older versions of PhysIO, a lot of questions have been discussed on our now defunct mailinglist 
   [tapas@sympa.ethz.ch](https://sympa.ethz.ch/sympa/info/tapas), 
   but its [archive](https://sympa.ethz.ch/sympa/arc/tapas) is still searchable.
3. For new requests, we would like to ask you to submit them as 
   [issues](https://github.com/translationalneuromodeling/tapas/issues) on our GitHub page for TAPAS.