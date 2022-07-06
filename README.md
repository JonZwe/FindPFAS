FindPFΔS
========

FindPFΔS (FindPolyFluoroDeltas) is an open source Python based algorithm
which can be used to extract high-resolution MS/MS spectra (HR-MS/MS)
that comprise certain fragment mass differences of interest (MS2 raw
data files). Its intended to be used for non-target screening for per-
and polyfluoroalkyl substances (PFAS). FindPFΔS is provided both as
Python source code and as Windows Application (.exe) with a graphical
user interphase (GUI).\
The executable version can be downloaded from Zenodo (https://doi.org/10.5281/zenodo.6797354) and can be
used independently of any programming language on a Windows computer.
The GUI was programmed with PySimpleGUI and the full source code which
was converted into an executable is also given.

Publication:\
Zweigle, J., Bugsel, B., Zwiener, C. (2022). "FindPFΔS: Non-target
screening for PFAS -- Comprehensive data mining for MS2 fragment mass
differences" submitted to Anal Chem

For further and more detailed information on functionality of FindPFΔS,
graphical user interface (GUI), output, and general recommendations also
refer to the publication and its Supporting Information.\
In case of questions regarding FindPFΔS: jozweigle\@gmail.com

General information:
--------------------

FindPFΔS was written in Python 3.9.7 and requires several Python
packages that can be found in the main script. In the following, the raw
data input and parameters that need to be specified for data evaluation
are briefly discussed. After explanation of each parameter, the
respective variables of the FindPFΔS source code are given in brackets.

### Parameters to search for MS2 fragment mass differences

-   Raw data input for FindPFΔS are MS2-files (.ms2). They can be
    generated vendor independently from most mass spectrometric data
    files by using the open-source software 'MSConvert' from the
    ProteoWizard toolbox (Download:
    https://proteowizard.sourceforge.io/index.html)

-   Prerequisites on raw data:

    1)  FindPFΔS only works with centroid data. If data was acquired in
        the profile mode, centroid data can be generated from profile
        data with the 'peak picking' function of MSConvert.
    2)  It is recommended to use FindPFΔS with only one collision energy
        (or e.g. one linear formula). If multiple collision energies
        were acquired each collision energy can be separatly converted
        to one MS2-file by using the 'subset' function of MSConvert.\
        (full\_path =
        'TestSample\_PFAS\_Standard\_MIX\_ddMS2\_20eV\_Inj5.ms2')

-   After choosing a datafile of interest, the fragment differences can
    be specified (e.g. CF2). Several differences can be specified at the
    same time. Note: Besides chemical formulae also exact masses can be
    used. In the GUI the input values have to be space separated. It is
    recommended to start by using only one difference per run to keep
    the results simple (and avoid to high false positive rates) in the
    beginning and to later use multiple differences.\
    (frags = \['CF2', 'C2F4', 'HF'\])

-   In the next step the number of differences desired for positive
    assignment can be specified: e.g. a difference has to be present at
    least 3 times per spectra (1 by default).\
    (n\_dists = 1)

-   Absolute mass tolerance (Da): It is recommended to start with a
    small tolerance (≤0.001 Da, depending on mass accuracy of the
    instrument) to keep the high false positive rates low in the
    beginning.\
    (tol = 0.001)

-   Next, the intensity threshold for fragments that should be
    considered can be choosen: (absolute or relative, in the GUI
    specified by Radio button). If the noise threshold of typical MS/MS
    spectra is known, absolute thresholds can have advantages, otherwise
    relative thresholds should be choosen.\
    (I\_min = 5)\
    (rel\_intens = True)

-   In the next field the mass tolerance (Da) to remove persistent
    background precursor masses that are triggered n-times for MS/MS can
    be set (see also next point).\
    (tol\_multiples = 0.002)

-   Number n above which precursor masses are excluded from data
    evaluation (e.g. if a compound is triggered 20 times for MS/MS it is
    likely a background signal without a peak shape) Note: The number 20
    is choosen arbitrarily. It strongly depends on the method parameters
    and the chromatograpic method. In cases where the focus is on high
    quality spectra and the instrument algorithm collects many spectra
    over a broad peak the number should be set rather high in the first
    place. If only e.g. \< 3 spectra over a peak are measured a number
    like e.g. 5 can be selected.\
    (n\_occur = 20)

-   Decision whether the two above mentioned criteria should be applied.
    If true (checked in GUI), multiple spectra from identical precursor
    masses are collected and only the spectrum with the highest
    precursor intensity will used for further evaluation (e.g. if a
    chromatographic peak has 5 MS/MS scans only the spectrum with the
    highest precursor intensity is extracted) Note: If unchecked, data
    evaluation will take much more time because every single
    MS2-spectrum in the dataset will be searched for fragment
    differences.\
    (match\_mzs = True)

-   Specification of sorting criteria (Radio button in GUI). It can be
    specified how the output data should be sorted (either by precursor
    intensity or number of fragments found positive). Both figures of
    spectra as well as the Excel output table
    'Results\_differences.xlsx' are sorted.\
    (sort\_intens = True)\
    (sort\_number = False)

-   Specifiy how many of the spectra with mass differences should be
    plotted.\
    (plot\_num = 20)

### Parameters for Kendrick mass defect analysis (homologous series)

-   Homologous series analysis will use the same repeating units as specified
    under fragment differences (e.g. CF2, CF2O or 115.9885)

-   Decision whether the MS1 masses should be screened for abundant
    repeating units (false (unchecked) by default; Note: If checked, the
    runtime can increase drastically).\
    (screen\_hs = False)

-   For MS1 peaks a Kendrick mass defect analysis will be performed.
    Mass tolerance to screen for homologous series.\
    (hs\_tol = 0.005)

-   Minimum number of homologues to be assigned as a series.\
    (n\_min = 3)
    
-   The output file is a KMD vs. m/z plot for visual inspection and 
    a [CF2]_homologous_series.xls file for in-depth analysis of homologous series.
    In the KMD plot, homologous series that have more homologues than specified
    in n_min are color coded based on the normalized retention time within their
    homologous series.
    The [CF2]_homologous_series.xls output file consists of a list with following entries:

	* Column A (number): numbering of features
	* m/z: measured m/z value by the instrument
	* RT: measured RT
	* mod: modulo, m/z mod [mass of repeating unit]
	* KMD: Kendrick mass defect
	* HS Number: Numbering of different homologous series. All the features with an identical HS Number vary by the mass of the repeating unit.
	* Homologoues: Number of features in the homologous series
	* Unique Homologues: Number of features in the homologous series, corrected for multiple entries of the same masses (e.g. a double integration of a split peak can lead to two features)
	* min Homologues: TRUE if Unique Homologues > n_min
	* Confirmed by MSMS: TRUE if also confirmed by fragment differences


### Optional: Parameters for suspect screening and diagnostic fragment screening

-   Optional file path for a suspect list to compare hits found by
    FindPFΔS with an underlying suspect list (the suspect list has to
    have the format used in suspect lists from the EPA CompTox Chemicals
    Dashboard
    (https://comptox.epa.gov/dashboard/chemical-lists/PFASOECD).
    Further, it has to be specified (Radio button in GUI) whether the
    data was acquired in positive or negative ionization mode. Note that
    either only \[M-H\]- or \[M+H\]+ adducts can be considered. Also
    the mass tolerance has to be defined.\
    (file\_susp\_list = False)\
    (tol\_suspects = 0.002)\
    (ESI\_neg = True)

-   Optional file for diagnostic fragments (.xlsx). All spectra detected
    by FindPFΔS will additionally be searched for the given diagnostic
    fragments.\
    (file\_diag\_list = False)

### Run FindPFΔS

-   When all mandatory parameters are specified, data evaluation with
    FindPFΔS can be started. It will generate a folder named after the
    sample. In this folder several Excel files with the results will be
    saved. Furthermore, a 'plots' folder is generated within the results
    folder with the generated figures.

Call for Contributions
----------------------

We appreciate help and suggestions from other reseachers to improve
FindPFΔS. Please don't hesitate to contact us.
