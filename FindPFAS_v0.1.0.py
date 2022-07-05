"""
FindPFAS (FindPolyFluoroDeltas) by Jonathan Zweigle (mass differences), Boris Bugsel (KMD)
Code to search MS2 data for fragment mass differences (DeltaFragments), e.g.:
(CF2)n | HF | C2F4+H2O | CF2O | CnH3F2n-3 | C2F4-HF = 79.9873 | C2F4-H2O = 81.9830 (diPAP, FTMAP, FTSA...)
and to find homologous series in MS1 data
"""
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from pylab import figure
import pyteomics
import os
import copy
from itertools import compress
from copy import deepcopy
from pyteomics import mass
from pyteomics.ms2 import IndexedMS2 # MS2
from pyteomics.mzml import MzML # mzML
from ismembertol import ismembertolerance # own function

# %% Data Input

# Input filepath to MS2-file here (.ms2)
full_path = 'TestSample_PFAS_Standard_MIX_ddMS2_20eV_Inj5.ms2'

# =============================================================================
# === PARAMETERS FOR MS2 FRAGMENT DIFFERENCE SEARCHING ========================
# =============================================================================
# Fragment & homologous series differences (Space separated)
frags = ['CF2', 'C2F4', 'HF']

# Minimum number of differences desired for positive assignment per spectrum
n_dists = 1

# Absolute mass tolerance in Dalton used for fragment mass difference searching
tol = 0.001

# Intensity threshold (absolute or relative) below which all MS2 peaks are removed
I_min = 5

# True for relative intensity
# False for absolute intensity
rel_intens = True

# Absolute mass tolerance in Dalton to search for identical precursor masses which are removed if they occur n times (n_occur)
tol_multiples = 0.002

# Number (n) of occurences above which identical precursor masses are removed before data evaluation
n_occur = 20

# Set to false only if all individual MS2-spectra should be searched. Note: Data evaluation will probably take much longer. If checked, only the spectrum with the highest precursor TIC from idetical masses is used for data evaluation
match_mzs = True

# Specify sorting criteria (by precursor intensity or number of fragments found with differences)
sort_intens = True
sort_number = False

# Number of plots desired
plot_num = 20

# =============================================================================
# === PARAMETERS FOR KENDRICK MASS DEFECT ANALYSIS (HOMOLOGOUS SERIES) ========
# =============================================================================
# Mass tolerance for homologous series (Da) 
hs_tol = 0.005
# Minimum homologues for valid homologous series
n_min = 3

# Check for automatic detection of abundant repeating units in MS data. Increases run time!
screen_hs = False


# =============================================================================
# === OPTIONAL PARAMETERS FOR SUSPECT SCREENING AND DIAGNOSTIC FRAGEMENTS =====
# =============================================================================

# OPTIONAL: input of suspect list (.xlsx) Note the suspect list format has to follow the structure of the EPA CompTox Chemicals Dashboard (https://comptox.epa.gov/dashboard/chemical-lists/PFASOECD)
file_susp_list = 'Chemical List PFASOECD-2022-06-14.xlsx' # 'Chemical List PFASOECD-2022-06-14.xlsx'
# Mass tolerance (Da) for suspect comparison
tol_suspects = 0.002
# Set to False if sample was measures in positive mode (ESI+)
ESI_neg = True

# =============================================================================
# OPTIONAL: input list with diagnostic fragments (.xlsx)
file_diag_list = False # 'Diagnostic_fragments_PFAS.xlsx'


# =============================================================================
# =============================================================================
# =============================================================================

# %% Convert fragments specified as numbers to float 
print('======================================\nFindPFAS started data evaluation...\n...')

for n in range(len(frags)):
    if frags[n][0].isnumeric():
        frags[n] = float(frags[n])


# %% Folder generation
# =========================================================================
# create results folder with sample name attached using os packages
sample_name = os.path.splitext(os.path.basename(full_path))[0]

Res_folder = 'Results_FindPFAS_' + sample_name
Plots_folder = os.path.join(Res_folder, 'plots')

# check if folder exists, if yes overwrite it
if not os.path.exists(Res_folder):
    os.makedirs(Res_folder)
if not os.path.exists(Plots_folder):
    os.makedirs(Plots_folder)
    
# =========================================================================
# create mass of difference array from chemical formulae
m_frags = np.array(np.zeros((len(frags))))
for n in range(len(frags)):
    if isinstance(frags[n], str):
        m_frags[n] = mass.calculate_mass(formula = frags[n])
    else:
        m_frags[n] = frags[n]

# %% Raw data reading
# =========================================================================
# check for file extension
[_, extension] = os.path.splitext(full_path)

if extension == '.ms2':
    print('Reading in MS2 data...\n...')
    
    def read_MS2(full_path):
        data = pyteomics.ms2.IndexedMS2(full_path) # read .ms2 file
        
        # loop to fill in m/z, RT, and intensity (TIC) of all precursor masses
        mz_list = []
        RT_list = []
        TIC_list = []
        spec_mz_list = []
        spec_intens_list = []
        for n in range(len(data)): 
            try:
                mz = data[n]['params']['precursor m/z']
            except KeyError:
                print('missing mz ignored ' + str(n))
                mz = False
            try:
                RT = data[n]['params']['RTime']
            except KeyError:
                print('missing RT ignored ' + str(n))
                RT = False
            try: 
                TIC = float(data[n]['params']['TIC']) # necessary to convert from str to float
            except KeyError:
                print('missing TIC ignored ' + str(n))
                TIC = False

            if mz and RT and TIC != False:
                mz_list.append(mz)
                RT_list.append(RT)
                TIC_list.append(TIC)
                spec_mz = data[n]['m/z array']
                spec_mz_list.append(spec_mz)
                spec_intens = data[n]['intensity array']
                spec_intens_list.append(spec_intens)
        
        mz_vec = np.array(mz_list) # convert lists to numpy arrays
        RT_vec = np.array(RT_list)
        TIC_vec = np.array(TIC_list)
        
        return mz_vec, RT_vec, TIC_vec, spec_mz_list, spec_intens_list
    
    mz_vec, RT_vec, TIC_vec, spec_mz_list, spec_intens_list = read_MS2(full_path)
    
    print(str(len(mz_vec)) + ' single MS2 spectra read...\n...')
#==========================================================================
elif extension == '.mzML':
    print('Reading in mzML data...\n...')
    
    def read_mzML(full_path):
        data = pyteomics.mzml.MzML(full_path)
        
        mz_list = []
        TIC_list = []
        RT_list = []
        spec_mz_list = []
        spec_intens_list = []
        # loop to read out each MS2 spectrum
        for n in range(len(data)):
            if data[n]['ms level'] == 2:
                mz_list.append(data[n]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
                TIC_list.append(data[n]['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity'])
                RT_list.append(data[n]['scanList']['scan'][0]['scan start time'])
                spec_mz_list.append(data[n]['m/z array'])
                spec_intens_list.append(data[n]['intensity array'])
    
        mz_vec = np.array(mz_list) # convert lists to numpy arrays
        RT_vec = np.array(RT_list)
        TIC_vec = np.array(TIC_list)
        
        return mz_vec, RT_vec, TIC_vec, spec_mz_list, spec_intens_list
    
    mz_vec, RT_vec, TIC_vec, spec_mz_list, spec_intens_list = read_mzML(full_path)
# =========================================================================

# %% FindPFAS Calculations
# =========================================================================
print("Remove persistent background precursor masses...\n...")
# loop to count the occurence of each m/z with a certain tolerance tol_multiples
count_vec = np.zeros(len(mz_vec))
for n in range(len(mz_vec)):
    counter = 0
    for m in range(len(mz_vec)):
        diff = abs(mz_vec[n] - mz_vec[m]) 
        if diff <= tol_multiples:
            counter = counter + 1      
    count_vec[n] = counter

unique_bool = count_vec <= n_occur # boolean with keep indices of unique m/z's

mz_vec_corr = mz_vec[unique_bool] # delete persistent background m/z's
RT_vec_corr = RT_vec[unique_bool]
TIC_vec_corr = TIC_vec[unique_bool]

# delete spectra of persistent background m/z's
spec_mz_list_corr = list(compress(spec_mz_list, unique_bool))
spec_intens_list_corr = list(compress(spec_intens_list, unique_bool))

# =========================================================================
# normalize spectra to basebeak if relative intensity threshold is checked
if rel_intens == True:
    for n in range(len(spec_intens_list_corr)):
        spec_intens_list_corr[n] = spec_intens_list_corr[n]/max(spec_intens_list_corr[n])*100

print(str(len(mz_vec_corr)) + ' single spectra used...\n...')

# =========================================================================
print("Remove intensities below threshold...\n...")
# intensity filtering
spec_intens_list_corr_fil = [None]*len(spec_intens_list_corr)
spec_mz_list_corr_fil = [None]*len(spec_mz_list_corr)
for n in range(len(spec_intens_list_corr)):
    keep_intens = spec_intens_list_corr[n] > I_min
    spec_intens_list_corr_fil[n] = spec_intens_list_corr[n][keep_intens]
    spec_mz_list_corr_fil[n] = spec_mz_list_corr[n][keep_intens]

# =========================================================================
# find unique m/z's and their maximum TIC and the corresponding spectrum
print("Match identical precursor masses...\n...")

if match_mzs == True:
    # find unique m/z's
    [uni_mzs, uni_idx] = np.unique(mz_vec_corr, return_index = True)

    # write the maximum TIC at the respective position of the uni_mzs array
    TIC_max = np.array(np.zeros((len(uni_mzs))))
    idx_TIC_max = [None]*len(uni_mzs)
    for n in range(len(uni_mzs)):
        TICs_ident_mz = []
        idx_TIC = []
        for m in range(len(mz_vec_corr)):
            if abs(uni_mzs[n] - mz_vec_corr[m]) == 0:
                TICs_ident_mz.append(TIC_vec_corr[m])
                idx_TIC.append(m)
        TIC_max[n] = max(TICs_ident_mz)
        idx_TIC_max[n]= idx_TIC[TICs_ident_mz.index(max(TICs_ident_mz))]

    # find RT and spectra at the respective positions
    RT_uni = RT_vec_corr[idx_TIC_max]
    
    # indexing list with a list of indices by mapping
    accessed_mapping_1 = map(spec_mz_list_corr_fil.__getitem__, idx_TIC_max)
    spec_mz_list_uni = list(accessed_mapping_1)
    
    accessed_mapping_2 = map(spec_intens_list_corr_fil.__getitem__, idx_TIC_max)
    spec_intens_list_uni = list(accessed_mapping_2)

    # rename arrays to original names
    mz_vec_corr = copy.deepcopy(uni_mzs)
    TIC_vec_corr = copy.deepcopy(TIC_max)
    RT_vec_corr = copy.deepcopy(RT_uni)
    spec_mz_list_corr_fil = copy.deepcopy(spec_mz_list_uni)
    spec_intens_list_corr_fil = copy.deepcopy(spec_intens_list_uni)

print(str(len(mz_vec_corr)) + ' spectra used...\n...')

# =========== Code to find MS/MS fragments with certain differences =======
# =========================================================================
# lower and upper boundary for difference between fragments
l_bound = m_frags - tol
u_bound = m_frags + tol

print('Calculating fragment mass differences...\n...')

def calc_dist_matrix(spec_mz_list_corr_fil):
    # create a list of matrices with n x n dimensions of each intensity array
    dist = [None]*len(spec_mz_list_corr_fil)
    for n in range(len(spec_mz_list_corr_fil)):
        dist[n] = np.array(np.zeros((len(spec_mz_list_corr_fil[n]),len(spec_mz_list_corr_fil[n]))))
    
    # fill list with the distances between fragments (n x n matrix)
    for s in range(len(spec_mz_list_corr_fil)):
        XX, YY = np.meshgrid(spec_mz_list_corr_fil[s], spec_mz_list_corr_fil[s])
        dist[s] = abs(XX-YY)
    
    return dist

dist = calc_dist_matrix(spec_mz_list_corr_fil)

# =========================================================================
# =========================================================================
print('Searching for fragment mass differences...\n...')

# create empty list for logicals
dist_log = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]

# loop which sets every value which is in the desired range to 1
# all other values are set to 0
for n in range(len(m_frags)):
    for m in range(len(spec_mz_list_corr_fil)):
        dist_log[n][m] = np.logical_and(dist[m] >= l_bound[n], dist[m] <= u_bound[n]).astype(int)                    


# create empty list
v = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]

# fill list with dimensions of each spectrum
for n in range(len(m_frags)): 
    for m in range(len(spec_intens_list_corr_fil)):
        v[n][m] = np.array(np.zeros(((len(spec_intens_list_corr_fil[m])))))
        
        
# loop that sums up all rows in each list (axis = 0)
# v is a boolean array which sets all values > 0 to 1
# values of 1 are the positions of the fragments which make up a distance
for n in range(len(m_frags)):
    for m in range(len(dist_log[0])):
        v[n][m] = np.array(np.sum(dist_log[n][m], axis = 0))
        v[n][m] = v[n][m] > 0        

#==========================================================================
# create empty lists: len(m_frags) x len(spectra_list)
spec_mz_pos = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]
spec_intens_pos = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]

# fill in zeros with the length of the positive peaks at the respective positions (sum(booleans))
for n in range(len(m_frags)):
    for m in range(len(spec_mz_list_corr_fil)):
        spec_mz_pos[n][m] =  np.array(np.zeros((sum(v[n][m]))))
        spec_intens_pos[n][m] = np.array(np.zeros((sum(v[n][m]))))
        
# list for each fragment mass difference with arrays with the positive fragments
# read out the spectra and save positive m/z's and intensities in new lists
for n in range(len(m_frags)):
    for m in range(len(spec_mz_list_corr_fil)):
        spec_mz_pos[n][m] = spec_mz_list_corr_fil[m][v[n][m]]
        spec_intens_pos[n][m] = spec_intens_list_corr_fil[m][v[n][m]]

# create list with number of fragments for each precursor m/z
n_frag = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]
for n in range(len(m_frags)):
    for m in range(len(spec_mz_list_corr_fil)):
        n_frag[n][m] = len(spec_mz_pos[n][m])


# sum of weighted intensities for all fragment differences
# NOT IN USE YET
sum_int_weighted = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_frags))]
for n in range(len(m_frags)):
    for m in range(len(spec_mz_list_corr_fil)):
        sum_int_weighted[n][m] = sum(spec_intens_pos[n][m]/len(spec_intens_pos[n][m]))


# sum up the total number of fragments for each precursor m/z
n_frag_tot = np.array(np.zeros((len(spec_mz_list_corr_fil))))
counter = 0
for n in range(len(spec_mz_list_corr_fil)):
    for m in range(len(n_frag)):
        counter = counter + n_frag[m][n]
        n_frag_tot[n] = counter
    counter = 0

idx_prec = n_frag_tot >= n_dists  # find all precursors which have n positive fragments
prec_total = sum(idx_prec) # number of hits

print(str(prec_total) + ' of ' + str(len(mz_vec_corr)) + \
      " m/z's were detected positive for fragment differences...\n...")

# =========================================================================
# ========== Saving data new for positive precursors ======================
# ========================================================================= 

# save positive precursors in new arrays
mz_pos = mz_vec_corr[idx_prec]
RT_pos = RT_vec_corr[idx_prec]
TIC_pos = TIC_vec_corr[idx_prec]

# save total number of positive fragments
n_frag_tot_pos = n_frag_tot[idx_prec]

# save new list with all precursors and their number of positive fragments
n_frag_pos = [[None] for x in range(len(m_frags))]
for n in range(len(m_frags)):
    n_frag_pos[n] = list(compress(n_frag[n], idx_prec))

# peaks with distances: delete all precursor with no positive fragments
spec_mz_POS = [[None] for x in range(len(m_frags))]
spec_intens_POS = [[None] for x in range(len(m_frags))]
for n in range(len(m_frags)):
    spec_mz_POS[n] = list(compress(spec_mz_pos[n], idx_prec))
    spec_intens_POS[n] = list(compress(spec_intens_pos[n], idx_prec))

# all peaks: delete all precur2sor with no positive fragments (new shorter variable name)
F1_pos = list(compress(spec_mz_list_corr_fil, idx_prec))
I1_pos = list(compress(spec_intens_list_corr_fil, idx_prec))

# %% Generate text file with parameters
# =========================================================================
# generate a text file with set parameters in Results folder
with open(os.path.join(Res_folder, 'FindPFAS_Parameters.txt'), "w") as file:
    file.write('Sample: ' + str(sample_name) + '\n' +
	       'Fragment differences: ' + str(frags)  + '\n' +
	       'Number of differences desired: ' + str(n_dists) + '\n' +
	       'Mass tolerance: ' + str(tol) + '\n' +
	       'Intensity threshold: ' + str(I_min) + '\n' +
	       'Remove multiples mass tolerance : ' + str(tol_multiples) + '\n' +
	       'Occurrence number threshold : ' + str(n_occur))

# %% Optional: Diagnostic fragments
# =========================================================================
# ===== OPTIONAL: Search for diagnostic fragments in positive spectra =====

if file_diag_list != False:
    
    print('Searching for diagnostic fragments...\n...')
    # read in excel file with diagnostic fragments
    dia_frags = pd.read_excel(file_diag_list)
    
    # keep only unique dia_frags and convert them to an array
    dia_frags_uni = np.array(dia_frags['F_containing_fragments_mass'])
    
    # create zeros list
    bool_dia_frags = [None]*len(F1_pos) 
    for n in range(len(bool_dia_frags)):
        bool_dia_frags[n] = np.array(np.zeros((len(F1_pos[n]))))

    # find diagnostic fragments in spectra
    for n in range(len(F1_pos)):
        [bool_dia_frags[n], _] = ismembertolerance(F1_pos[n], dia_frags_uni, tol)

    # create empty list with zeros for positive dia_frags
    spec_mz_dias = [None]*len(F1_pos) # create zeros list
    for n in range(len(spec_mz_dias)):
        spec_mz_dias[n] = np.array(np.zeros((sum(bool_dia_frags[n]))))
    
    # copy list for intensities
    spec_intens_dias = copy.deepcopy(spec_mz_dias)
    
    # save positive diagnostic fragments in list with arrays
    for n in range(len(spec_mz_dias)):
        spec_mz_dias[n] = F1_pos[n][bool_dia_frags[n]]
        spec_intens_dias[n] = I1_pos[n][bool_dia_frags[n]]
    
    # count diagnostic fragments
    n_dia_frags = np.array(np.zeros((len(spec_mz_dias))))
    for n in range(len(n_dia_frags)):
        n_dia_frags[n] = sum(bool_dia_frags[n])
    
    # find indizes of m/z's with diagnostic fragments
    idx_dia_frags = n_dia_frags > 0
    n_dia_frag_prec_tot = sum(idx_dia_frags)
    
    print(str(n_dia_frag_prec_tot) + ' of ' + str(len(mz_pos)) + \
          " m/z's have diagnostic fragments...\n...")
    
    # save precursors with diagnostic fragments in new arrays
    mz_pos_dias = mz_pos[idx_dia_frags]
    RT_pos_dias = RT_pos[idx_dia_frags]
    TIC_pos_dias = TIC_pos[idx_dia_frags]
    n_dia_frag_prec = n_dia_frags[idx_dia_frags]
    
    # save spectra with diagnostic fragments in new arrays
    spec_mz_dias_pos = list(compress(spec_mz_dias, idx_dia_frags))
    spec_intens_dias_pos =  list(compress(spec_intens_dias, idx_dia_frags))
    
    # all spectra peaks
    F1_pos_dias = list(compress(F1_pos, idx_dia_frags))
    I1_pos_dias = list(compress(I1_pos, idx_dia_frags))
    
    # create pandas DataFrame for m/z's with diagnostic fragments
    D_dias = {'m/z': mz_pos_dias,'RT': RT_pos_dias,'TIC': TIC_pos_dias, 'n_tot': n_dia_frag_prec, \
              'mz_peaks': F1_pos_dias, 'intens_peaks': I1_pos_dias, \
              'mz_peaks_diagnostic': spec_mz_dias_pos, 'intens_peaks_diagnostic': spec_intens_dias_pos}
    Df_dias = pd.DataFrame(data = D_dias)
    
    # sort DataFrame according to number of diagnostic fragments
    Df_s_dias = Df_dias.sort_values(by='n_tot',ascending = False)
    Df_s_dias = Df_s_dias.reset_index(drop = True)

    # write DataFrame from diagnostic fragments to Excel file
    Df_s_dias.to_excel(os.path.join(Res_folder, "Results_diagnostic_fragments.xlsx"))
    
# %% Optional Suspect screening
# =========================================================================
# ======== OPTIONAL: Suspect list comparison ==============================
if file_susp_list != False:
    
    print("Compare positive precursor m/z values with suspect list...\n...")
    # read in suspect list from Excel file
    susp_list = pd.read_excel(file_susp_list).dropna(subset=['AVERAGE MASS'])
    mono_mass = np.array(susp_list['MONOISOTOPIC MASS'].astype(dtype ='float64')) # monoisotopic mass array
    mass_H = 1.007276 # mass of hydrogen
    if ESI_neg == True:
        mz_pos_ESI_H = mz_pos + mass_H # add hydrogen -> convert to mass (only [M-H]-)
    else:
        mz_pos_ESI_H = mz_pos - mass_H # substract hydrogen -> convert to mass (only [M+H]+)
    
    # compare masses of hits with monoisotopic masses in suspect list
    [hit_mz_sus,hit_susp] = ismembertolerance(mz_pos_ESI_H, mono_mass, tol_suspects)
    
    # read out positive hits
    mz_suspect_hit = mz_pos[hit_mz_sus]
    RT_suspect_hit = RT_pos[hit_mz_sus]
    TIC_suspect_hit = TIC_pos[hit_mz_sus]
    
    CpdName_hit = susp_list['PREFERRED NAME'][hit_susp] # compound name
    CAS_hit = susp_list['CASRN'][hit_susp] # CAS number
    Mono_mass_hit = susp_list['MONOISOTOPIC MASS'][hit_susp] # monoisotopic mass of hits (from suspect list)
    
    print(str(len(CpdName_hit)) + ' of ' + str(prec_total) + \
          " m/z values have an accurate mass match with the suspect list...\n...")
    
    # create pandas DataFrame with hits from suspect list
    D_susp_exp = {'m/z': mz_suspect_hit,'RT': RT_suspect_hit,'TIC': TIC_suspect_hit}
    D_susp = {'mono_mass': Mono_mass_hit,'CAS': CAS_hit,'CpdName': CpdName_hit}
    Df_susp_exp = pd.DataFrame(data = D_susp_exp)
    Df_susp = pd.DataFrame(data = D_susp)
    
    # write DataFrame of suspect hits to Excel file
    Df_susp.to_excel(os.path.join(Res_folder, "Results_suspect.xlsx"))
    Df_susp_exp.to_excel(os.path.join(Res_folder, "Results_suspect_experimental.xlsx"))

# %% Save final data
# =========================================================================
# ========= CREATE ONE DATAFRAME FOR FRAGMENT MASS DIFFERENCE HITS ========

# create pandas DataFrame for general data
D = {'m/z': mz_pos,'RT': RT_pos,'TIC': TIC_pos, 'n_tot': n_frag_tot_pos, \
      'mz_peaks': F1_pos, 'intens_peaks': I1_pos}
Df_general = pd.DataFrame(data = D)


# convert fragment masses to string for column names
frags_str = [[None] for x in range(len(frags))]
for n in range(len(frags)):
    frags_str[n] = str(frags[n])

# create pandas DataFrame for n_frag_pos
Df_frags = pd.DataFrame(n_frag_pos).T # transpose
Df_frags.columns = frags_str
# IDEA: MAKE PCA OF Df_pos_frags TO FIND OUT WHICH FRAGMENTS RELATED TO EACH OTHER

# add prefix to fragment string
append_str_frag = 'mz_'
append_str_intens = 'intens_'
spec_frag_str = [append_str_frag + sub for sub in frags_str]
spec_intens_str = [append_str_intens + sub for sub in frags_str]

# create pandas DataFrames for positive spectra peaks
Df_pos_spec_mz = pd.DataFrame(spec_mz_POS).T
Df_pos_spec_mz.columns = spec_frag_str
Df_pos_spec_intens = pd.DataFrame(spec_intens_POS).T
Df_pos_spec_intens.columns = spec_intens_str

# merge DataFrames together
Df = pd.concat([Df_general, Df_frags.reindex(Df_general.index),\
                Df_pos_spec_mz.reindex(Df_general.index),\
                Df_pos_spec_intens.reindex(Df_general.index)],\
                axis = 1)

# condition for sorting final reuslts
if sort_intens == True:
    Df_s = Df.sort_values(by='TIC',ascending=False)
    Df_s = Df_s.reset_index(drop=True)
else:
    Df_s = Df.sort_values(by='n_tot',ascending=False)
    Df_s = Df_s.reset_index(drop=True)
    
# write DataFrame from differences to Excel file
Df_s.to_excel(os.path.join(Res_folder, "Results_differences.xlsx"))


# %% Plotting
#==========================================================================
#========================  PLOTTING RESULTS ===============================
# plotting fragments with mass differences

colors = ['deeppink','mediumpurple','indianred','orangered', 'seagreen','deepskyblue','royalblue','limegreen','orange','cyan','gold','blue']
# NOTE: max 12 colors -> max 12 differences above -> ERROR!
count = 0
for s in range(len(frags_str)):
    # sort DataFrame according to fragment input
    Df_p = None
    Df_p = Df_s[Df_s[frags_str[s]] != 0]
    Df_p = Df_p.reset_index(drop=True)
    
    if Df_p.empty == True:
        n_plots = 0
    
    else:
        if sort_intens == True:
            Df_p = Df_p.sort_values('TIC', ascending = False, ignore_index = True)
        else:
            Df_p = Df_p.sort_values(frags_str[s], ascending = False, ignore_index = True)
        
        if len(Df_p[[frags_str[s]]]) < plot_num and len(Df_p[[frags_str[s]]]) > 0:
            n_plots = len(Df_p[[frags_str[s]]])
        else:
            n_plots = plot_num


    for n in range(n_plots):
        count += 1
        fig1 = figure(count,figsize=(12,6))
        plt.stem(Df_p['mz_peaks'][n],Df_p['intens_peaks'][n],'k',\
                  markerfmt=" ",basefmt=" ")
    
        plt.stem(Df_p['mz_' + frags_str[s]][n],Df_p['intens_' + frags_str[s]][n],colors[s],\
                  markerfmt=" ",basefmt=" ")
            
        for i, txt in enumerate(Df_p['mz_' + frags_str[s]][n]):
            plt.annotate(txt, (Df_p['mz_' + frags_str[s]][n][i],Df_p['intens_' + frags_str[s]][n][i]), color = colors[s], rotation = 20, fontsize=7)
            
        markerline, stemlines, baseline = \
            plt.stem(Df_p['mz_' + frags_str[s]][n],\
                      Df_p['intens_' + frags_str[s]][n],colors[s],markerfmt=" ",basefmt=" ")
        plt.setp(stemlines, 'linewidth', 3)
        plt.ylim(ymin = 0)
        
        plt.title('\u0394' +  frags_str[s] + ', m/z = ' + str(Df_p['m/z'][n]) + \
                  ', RT = ' + str(Df_p['RT'][n]) + \
                  ', TIC = ' + str(Df_p['TIC'][n]) + \
                  ', n_frags = ' + str(Df_p[frags_str[s]][n]), fontsize=10)
        plt.xlabel('m/z',size=20)
        plt.ylabel('Counts',size=20)
        plt.xticks(fontsize=18) # change font size of xticks
        plt.yticks(fontsize=18)

        fig_name = frags_str[s] + '_spectrum_' + str(n+1) + '.png'
        fig1.savefig(os.path.join(Plots_folder, fig_name), bbox_inches='tight',\
                    format='png', dpi=100)
        plt.close(fig1)
        
# =========================================================================  
# =========================================================================
if file_diag_list != False: # only plot if diagnostic fragments are there
# plotting diagnostic fragments
    n_dias = 20
    if len(Df_s_dias) < n_dias:
        n_dias = len(Df_s_dias)
    else:
        n_dias = 20
        
    for n in range(n_dias):
        fig2 = figure(n,figsize=(12,6))
        plt.stem(Df_s_dias['mz_peaks'][n],Df_s_dias['intens_peaks'][n],'k',\
                  markerfmt=" ",basefmt=" ")
    
        plt.stem(Df_s_dias['mz_peaks_diagnostic'][n],Df_s_dias['intens_peaks_diagnostic'][n],'b',\
                  markerfmt=" ",basefmt=" ")
            
        for i, txt in enumerate(Df_s_dias['mz_peaks_diagnostic'][n]):
            plt.annotate(txt, (Df_s_dias['mz_peaks_diagnostic'][n][i],Df_s_dias['intens_peaks_diagnostic'][n][i]), color = 'blue', rotation = -30, fontsize=8)
            
            
        markerline, stemlines, baseline = \
            plt.stem(Df_s_dias['mz_peaks_diagnostic'][n],Df_s_dias['intens_peaks_diagnostic'][n],'b',\
                      markerfmt=" ",basefmt=" ")
        plt.setp(stemlines, 'linewidth', 3) 
        
        plt.title('m/z = ' + str(Df_s_dias['m/z'][n]) + ', RT = ' + str(Df_dias['RT'][n]) + \
                  ' TIC = ' + str(Df_s_dias['TIC'][n]))
        plt.xlabel('m/z',size=20)
        plt.ylabel('Counts',size=20)
        plt.xticks(fontsize=18) # change font size of xticks
        plt.yticks(fontsize=18)
        plt.ylim(ymin = 0)
        
        fig_name = 'spectrum_diagnostic' + str(n+1) + '.png'
        fig2.savefig(os.path.join(Plots_folder, fig_name), bbox_inches='tight',\
                    format='png', dpi=100)
        plt.close(fig2)

# =========================================================================
# Heatmap with number of fragment differences
fig3, ax = plt.subplots()
im = ax.imshow(Df_frags.T, cmap='viridis',aspect='auto')
ax.set_xticks(np.arange(len(Df_s['m/z'])))
ax.set_xticklabels(labels=Df_s['m/z'])
ax.set_yticks(np.arange(len(frags_str)))
ax.set_yticklabels(labels = frags_str)
fig3.colorbar(im, ax=ax)

plt.setp(ax.get_xticklabels(), rotation = 90, ha="right",
          rotation_mode="anchor",size = 7)
ax.set_title('Number of fragment differences')
plt.xlabel('m/z',size=10)
plt.ylabel('Fragment difference',size=10) 
fig3.savefig(os.path.join(Plots_folder, 'heatmap.png'), bbox_inches='tight',\
            format='png', dpi=200)
plt.close(fig3)

# %% Calculation of frequent MS1 mass differences
# =========================================================================
# Create homologous series (HS) that vary by the same mass differences.
# For n triggered precursor compounds, an n-by-n matrix is generated with
# differences from any compound to any compound. The differences are then 
# sorted according to their occurrence. This procedure allows the automatic
# detection of abundant repeating units.
# =========================================================================

if screen_hs == True:   
    print('Find abundant repeating units... \n...')
    # Ceate zero array for all precursors
    dist_hs = np.array(np.zeros((len(mz_vec_corr),len(mz_vec_corr))))

    # fill dist_hs array with differences of all precursors
    NN, MM = np.meshgrid(mz_vec_corr, mz_vec_corr)
    dist_hs = abs(NN-MM)

    # Convert distances in matrix to linear array
    dist_hs_vec = np.concatenate(dist_hs, axis = 0)

    # Remove distances outside boundaries (default: 18 - 250)
    dist_hs_vec_corr = dist_hs_vec[np.logical_and(dist_hs_vec > 18, dist_hs_vec < 250)]

    # Sort dist_hs ascending and save to dist_hs_sorted
    dist_hs_sorted = sorted(dist_hs_vec_corr)

    # Check occurrence of each value and find abundant repeating units
    avg_mass_temp = []
    avg_mass = []
    count_hs = np.zeros(len(dist_hs_sorted))
    for n in range(len(dist_hs_sorted)-1):
        counter = 1
        i = n
        avg_mass_temp = [dist_hs_sorted[n]]
        while dist_hs_sorted[i+1]-dist_hs_sorted[n] < 0.01 and i < len(dist_hs_sorted)-2:
            counter = counter + 1
            i = i+1
            count_hs[n] = counter
            avg_mass_temp.append(dist_hs_sorted[i])
        avg_mass.append(avg_mass_temp)
        avg_mass_temp = []
    avg_mass.append([dist_hs_sorted[-1]])


    # Create del_index to delete entries from same HS
    del_index = []
    for n in range(len(count_hs)-1):
        if count_hs[n+1] == count_hs[n]-1 or count_hs[n+1] == 0:
            del_index.append(n+1)

    count_arr = np.delete(count_hs, del_index)
    dist_hs_vec_corrected = np.delete(dist_hs_sorted, del_index)
    avg_mass = np.asarray(avg_mass,dtype='object') ################## Error message????
    np.delete(avg_mass, del_index)

    # Create Dataframe with detected repeating units
    avg_mass_list = [i for j, i in enumerate(avg_mass) if j not in del_index] 
    avg_masses = []
    for n in range(len(avg_mass_list)):
        avg_masses.append(sum(avg_mass_list[n])/len(avg_mass_list[n]))

    M = {'mz':dist_hs_vec_corrected,'average mz':avg_masses,'counts':count_arr}
    HS_Dataframe = pd.DataFrame(data=M)
    HS_Dataframe.sort_values('counts', ascending=False, inplace=True, ignore_index=True) #sort according to abundance
    
    # Write output file
    Df_HS_red = HS_Dataframe.drop('mz', axis=1)
    Df_HS_red.to_excel(os.path.join(Res_folder, 'Abundant_repeating_units.xlsx'))

# %% Homologous series detection
# =========================================================================
# Create homologous series (HS) based on the user input
# =========================================================================

mz_HS_tot = np.array([]) # array with all m/z's found by HS for functions
for j in range(len(m_frags)):
    print('Create KMD plot for', m_frags[j], '... \n...')
    rep_unit = m_frags[j]
    
    #Calculate modulo: Compounds from the same homologous series bear an identical modulo
    modulo = mz_vec_corr % rep_unit

    # Calculation of Kendrick masses and features that are in the same homologous series
    KM = mz_vec_corr*round(rep_unit)/rep_unit
    KM_round = np.round_(KM)
    KMD = KM - KM_round
    Mod_HS = {'mz':mz_vec_corr,'RT':RT_vec_corr,'mod':modulo,'KMD':KMD}
    Mod_HS_Dataframe = pd.DataFrame(data=Mod_HS)
    Mod_HS_Dataframe.sort_values('mod', inplace=True, ignore_index=True) #nach aufsteigendem modulo sortieren

    HS_num = np.zeros(len(mz_vec_corr))
    # HS_abund_corr = np.zeros(len(mz_vec_unique))
    exit_loop = 0
    #HS_num[0] = 1
    for n in range(len(mz_vec_corr)-1):
        i = n
        if HS_num[n] == 0 and n != 0:
            HS_num[n] = n
            while Mod_HS_Dataframe['mod'][i+1]-Mod_HS_Dataframe['mod'][n] < hs_tol and exit_loop == 0: #and i < (len(mz_vec_unique)-2):
                HS_num[i+1] = n
                HS_num[n] = n
                if i < (len(mz_vec_corr)-2):
                    i = i+1
                else:
                    exit_loop = 1
    if HS_num[-1] == 0:
        HS_num[-1] = n+1

    Mod_HS_Dataframe['HS Number'] = HS_num

    # Calculate number of members in HS
    HS_num_temp = np.unique(HS_num, return_counts = True)
    hsnumber = []
    for n in range(len(HS_num_temp[1])):
        for s in range(HS_num_temp[1][n]):
            hsnumber.append(HS_num_temp[1][n])

    Mod_HS_Dataframe['Homologues'] = hsnumber
    
    # Calculate corrected number of members in HS (same masses are neglected)
    HS_num2 = np.zeros(len(mz_vec_corr))
    for n in range(len(HS_num2)):
        HS_num2[n] = len(np.unique(round(Mod_HS_Dataframe['mz'][Mod_HS_Dataframe['HS Number'] == Mod_HS_Dataframe['HS Number'][n]])))
        
    Mod_HS_Dataframe['Unique Homologues'] = HS_num2

    # Check if number of homologues are greater than specified value
    HS_min_check = []
    for n in range(len(hsnumber)):
        if Mod_HS_Dataframe['Unique Homologues'][n] >= n_min:
            HS_min_check.append(True)
        else:
            HS_min_check.append(False)

    Mod_HS_Dataframe['min Homologues'] = HS_min_check

    # Check the features that were also detected by their fragments
    Mod_HS_Dataframe.sort_index(inplace=True)

    MSMS_confirmed = []
    for n in range(len(Mod_HS_Dataframe['mz'])):
        check = Mod_HS_Dataframe['mz'][n] in mz_pos
        if check == True:
            MSMS_confirmed.append(True)
        elif check == False:
            MSMS_confirmed.append(False)

    # Attach to Mos_HS_Dataframe
    Mod_HS_Dataframe['Confirmed by MSMS'] = MSMS_confirmed

    # Create KMD plots and output files
    Mod_HS_Dataframe.sort_values('min Homologues', inplace=True) #nach aufsteigendem modulo sortieren
    colors = {True:'fuchsia', False:'gainsboro'}
    Mod_HS_Dataframe.sort_values('Confirmed by MSMS', inplace=True) #nach aufsteigendem modulo sortieren

    fig1 = figure(count,figsize=(8,6))

    # Plot grey scatters (feature not confirmed HS)
    Mod_HS_Dataframe_grey = Mod_HS_Dataframe[~(Mod_HS_Dataframe['min Homologues'])]
    plt.scatter(Mod_HS_Dataframe_grey['mz'], Mod_HS_Dataframe_grey['KMD'], c='gainsboro', label='Feature')

    # Plot colored scatters & calculate normalized RT: Confirmed by HS
    Mod_HS_Dataframe_pos = Mod_HS_Dataframe[(Mod_HS_Dataframe['min Homologues'])] #Mod_HS_Dataframe[(Mod_HS_Dataframe['min Homologues'] | Mod_HS_Dataframe['Confirmed by MSMS'])]

    HS_unique = np.unique(Mod_HS_Dataframe_pos['HS Number'])
    RT_norm_fin = pd.Series(dtype='float')
    for i in range(len(HS_unique)):
        HS_temp = Mod_HS_Dataframe_pos[Mod_HS_Dataframe_pos['HS Number'] == HS_unique[i]]
        RT_min = np.amin(HS_temp['RT'])
        RT_max = np.amax(HS_temp['RT'])
        RT_norm = ((HS_temp['RT'] - RT_min) / (RT_max - RT_min))
        RT_norm_fin = RT_norm_fin.append(RT_norm)
    Mod_HS_Dataframe_pos = pd.concat([Mod_HS_Dataframe_pos,RT_norm_fin.rename('RT_norm3')],axis=1)

    plt.scatter(x=Mod_HS_Dataframe_pos['mz'], y=Mod_HS_Dataframe_pos['KMD'], c=np.array(Mod_HS_Dataframe_pos['RT_norm3']), cmap='viridis', label='Homologous series')
    plt.colorbar(label='Normalized Retention Time')
    
    # Plot red circles: Confirmed by MSMS
    Mod_HS_Dataframe_MSMS = Mod_HS_Dataframe[(Mod_HS_Dataframe['Confirmed by MSMS'])]
    plt.scatter(Mod_HS_Dataframe_MSMS['mz'], Mod_HS_Dataframe_MSMS['KMD'], facecolors='none', edgecolors="red", label='MSMS hit', s=100, alpha=0.3)

    # Create labels
    #title = round(rep_unit, 4)
    plt.title('Repeating Unit ' +frags_str[j])
    plt.xlabel('Precursor m/z')
    plt.ylabel(frags_str[j]+'-based KMD')
    plt.legend()

    fig_name = frags_str[j] + '_homologous_series.png'
    fig1.savefig(os.path.join(Plots_folder, fig_name), bbox_inches='tight',\
                format='png', dpi=100)
    plt.close(fig1)
    
    # Write HS output file
    Mod_HS_Dataframe_pos_cp = Mod_HS_Dataframe_pos.copy()
    Mod_HS_Dataframe_pos_sort = Mod_HS_Dataframe_pos.sort_values(by='HS Number')
    Mod_HS_Dataframe_pos_sort.to_excel(os.path.join(Res_folder, frags_str[j] + '_homologous_series.xlsx'))
    
    # append all m/z's found by HS
    mz_HS_tot = np.append(mz_HS_tot,np.array(Mod_HS_Dataframe_pos_sort['mz']))
    
# save only unique m/z's    
mz_HS_tot = np.unique(mz_HS_tot)    
#=================================HSend======================================== 

print('Evaluation successfully finished!\n')
