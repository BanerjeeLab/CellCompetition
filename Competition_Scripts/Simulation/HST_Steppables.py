# v4.5 only needs one cc3d import
from cc3d.core.PySteppables import *

import numpy as np
import random as random
import csv as csv
import time as time
import os as os
import sys
import json as js
  

######################################################################
##             Fundamental Paramters of the Simulation              ##
######################################################################
global CELL_ATTR_TRACKER
CELL_ATTR_TRACKER = {}
CELL_ATTR_TRACKER['MITOSIS_COND_DICT'] =    {
                                                        'sizer_thresholds':{}, 
                                                        'timer_thresholds':{}, 
                                                        'timer_start_mcs' :{} 
                                                        } 
CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'] = {}
CELL_ATTR_TRACKER['COLOR_CELLS'] = {}      

SIM_CONSTS = {}

SIM_CONSTS['SIM'] = {
    'X_MAX': 1500,      ## Match value in xml file
    'Y_MAX': 1500,      ## Match value in xml file
    'EDGE_BUFFER': 50,  
    
    'MUM2_TO_PIXELS'  : 0.33**(-2.), # [um^2]  * 0.33^-2 [pix/um^2]  =  [pix]
    'PIX_TO_MUM2'     : 0.33**( 2.), # [pix]   * 0.33^2  [um^2/pix]  =  [um^2] 
    'MUM_TO_1D_LENGTH': 0.33,        # [1D SIM LENGTH] * 0.33 [um]/[1D SIM LENGTH] = [um]

    'SAVE_DATA_TO_CSV_INTERVAL': 500,   # Save data every  X  mcs
    'END_SIM_IF_NO_CELLS'      : True,
    
    'P_MECH_MIN_PER_FRAME': 4. ## See Bove et.al. 2017 "Local cellular neighborhood controls proliferation in cell competition"
    }
    

##########################################################
### Where to save the data (LINUX file system)
###
### NOTE: This path must point to the folder Compucell3D 
###       makes at the start of the simulation. Otherwise
###       permission errors can/will occur.
###
###       If permissions issues present a problem you can  
###       stop all data from saving with the edit:
###          SIM_CONSTS['SIM']['SAVE_DATA'] = False 
###
##########################################################
SIM_CONSTS['SIM']['DATA_OUT_PATH'] = '/home/export/lcarpent/CompuCell3D/output/' 
# Note WINDOWS file systems need
#SIM_CONSTS['SIM']['DATA_OUT_PATH'] = 'C:\\...\\...\\'

SIM_CONSTS['SIM']['SAVE_DATA'] = True
Ave_Vol_Data
Cell_Heritage_Data
Cell_Cycle_Data
Local_Density_Data
P_Apo_Density_Data
Contact_Death_Data
Vol_Den_Count_Data
Extrusion_Data
Tot_Div_Apo_Count_Data



## A couple optional conditions for ending a simulation
KILL_SIM_IF_MUT_ELIMINATED = True
KILL_SIM_IF_WT_ELIMINATED  = True
KILL_SIM_IF_MUT_EXCEEDS_WT = False

KILL_SIM_IF_R_MUT_EXCEEDS_THRESHOLD = False 
R_MUT_THRESHOLD = 700


## Choose which forms a cell death are calculated (left on by default)
SIM_CONSTS['SIM']['CALCULATE_MECH_PROB_APO']    = True
SIM_CONSTS['SIM']['CALCULATE_BIOCHEM_PROB_APO'] = True
    
######################################################################
##                General Parameters of the Sim                     ##
######################################################################    

## Decide which plot you want the cc3d to show
SIM_CONSTS['SIM']['SHOW_DENSITY_APO_COUNT_PLOT']  = False
SIM_CONSTS['SIM']['SHOW_CONTACT_APO_COUNT_PLOT']  = False
SIM_CONSTS['SIM']['SHOW_DIVISION_COUNT_PLOT']     = False
SIM_CONSTS['SIM']['SHOW_POPULATION_COUNT_PLOT']   = True
SIM_CONSTS['SIM']['SHOW_AVE_ADDED_dATdt_PLOT']    = True
SIM_CONSTS['SIM']['SHOW_AVE_AREA_PLOT']           = True
SIM_CONSTS['SIM']['SHOW_AVE_DENSITY_PLOT']        = True
SIM_CONSTS['SIM']['SHOW_WT_DIV_APO_TOTALS_PLOT']  = True
SIM_CONSTS['SIM']['SHOW_MUT_DIV_APO_TOTALS_PLOT'] = True

######################################################################
##                                                                  ##
######################################################################
SIM_CONSTS['WT']  = {}

# Choose a hardcoded SC ave area
SIM_CONSTS['WT']['subconfluent average area']     = 3200 # [pix]

SIM_CONSTS['WT']['simulation initial target area'] = 0.5 * SIM_CONSTS['WT']['subconfluent average area']
SIM_CONSTS['WT']['simulation initial cell length'] = 25 # [pix],  cell starts as a square

SIM_CONSTS['WT']['timer threshold'] = 225 # [mcs]
SIM_CONSTS['WT']['timer threshold noise'] = 0.025 * SIM_CONSTS['WT']['timer threshold']

SIM_CONSTS['WT']['sizer threshold']       = 1400 # pix
SIM_CONSTS['WT']['sizer threshold noise'] =    0.025 * SIM_CONSTS['WT']['sizer threshold']

SIM_CONSTS['WT']['contact inhibition'] = 0.01

###  To investigate the effect of cell stiffness being drawn from
#       a distribution. I never got around to it so I set std = 0.
SIM_CONSTS['WT']['Cell Stiffness Skew']        = -2.    # 
SIM_CONSTS['WT']['Cell Stiffness Mean']        =  1.    # 
SIM_CONSTS['WT']['Cell Stiffness std' ]        =  0.    # To set all cells stiffnesses equal to the mean set this to 0.
SIM_CONSTS['WT']['Cell Stiffness Lower Bound'] =  0.    # 
SIM_CONSTS['WT']['Cell Stiffness Upper Bound'] =  2.    # 

SIM_CONSTS['WT']['growth rate'] = 0.5 * SIM_CONSTS['WT']['subconfluent average area'] / SIM_CONSTS['WT']['timer threshold']   #  [pix/mcs]

SIM_CONSTS['WT']['density p_apo,max']            =   0.0015
SIM_CONSTS['WT']['density death sensitivity']    = 145.2     * SIM_CONSTS['SIM']['MUM2_TO_PIXELS']
SIM_CONSTS['WT']['density death susceptibility'] =   0.01001 * SIM_CONSTS['SIM']['PIX_TO_MUM2']

SIM_CONSTS['WT']['biochemical p_apo,max'] = 1.
SIM_CONSTS['WT']['biochemical steepness'] = 1.35
SIM_CONSTS['WT']['biochemical Hill coef'] = 3.6   

SIM_CONSTS['WT']['extrusion threshold'] = 0.25


######################################################################
## MUT defaults are a clone of WT parameters
SIM_CONSTS['MUT']  = {}

for k in SIM_CONSTS['WT']:
    SIM_CONSTS['MUT'][k] = SIM_CONSTS['WT'][k]


###################################
##  Targeted Phase Space Search  ##
###################################

Prob_Contact_Death = {  'P_bio_050': { 'p_max': 0.95 , 'Steepness': 1.65, 'Hill coef': 3.6  } ,
                        'P_bio_100': { 'p_max': 1.   , 'Steepness': 1.35, 'Hill coef': 3.6  } ,
                        'P_bio_125': { 'p_max': 0.975, 'Steepness': 1.45, 'Hill coef': 2.75 } ,
                        'P_bio_150': { 'p_max': 0.8  , 'Steepness': 1.45, 'Hill coef': 2.1  } ,
                        'P_bio_175': { 'p_max': 0.55 , 'Steepness': 1.1 , 'Hill coef': 1.85 } ,
                        'P_bio_200': { 'p_max': 0.4  , 'Steepness': 0.75, 'Hill coef': 1.75 } }

# Contact based competition prob_death {{Contact_Death_Prob}} set in ParameterScanSpecs.json
p_con_key = {{ Contact_Death_Prob }}
SIM_CONSTS['MUT']['biochemical p_apo,max'] = Prob_Contact_Death[p_con_key]['p_max']
SIM_CONSTS['MUT']['biochemical steepness'] = Prob_Contact_Death[p_con_key]['Steepness']
SIM_CONSTS['MUT']['biochemical Hill coef'] = Prob_Contact_Death[p_con_key]['Hill coef'] 

# mut_k_lam_vals = [  [0.0010221, 0.2474125] , [0.0009784, 0.4425875] , [0.0015470, 0.5450875] , [0.0016160, 0.3499125] , [0.0024460, 0.6475875] ,
                    # [0.0025552, 0.4524125] , [0.0008966, 0.8329374] , [0.0014176, 0.9354374] , [0.0022415, 1.0379374] , [0.0027884, 0.0620626] ,
                    # [0.0008216, 1.2232874] , [0.0012991, 1.3257874] , [0.0020540, 1.4282874] , [0.0007612, 1.5648436] , [0.0012035, 1.6673436] ,
                    # [0.0019029, 1.7698436] , [0.0010790, 2.1552810] , [0.0017061, 2.2577810] , [0.0009674, 2.6432184] , [0.0015297, 2.7457184] ,
                    # [0.0008674, 3.1311558] , [0.0013715, 3.2336558] , [0.0034389, 0.7460804] , [0.0036349, 0.5539196] , [0.0048633, 0.8460804] ,
                    # [0.0051406, 0.6539196] , [0.0030779, 1.1304019] , [0.0040612, 0.1695981] , [0.0043528, 1.2304019] , [0.0057434, 0.2695981] ,
                    # [0.0027548, 1.5147234] , [0.0038959, 1.6147234] , [0.0025001, 1.8510047] , [0.0035357, 1.9510047] , [0.0021765, 2.3314066] ,
                    # [0.0030780, 2.4314066] , [0.0018948, 2.8118085] , [0.0026796, 2.9118085] , [0.0016495, 3.2922104] , [0.0023328, 3.3922104] ,
                    # [0.0072442, 0.9937878] , [0.0077648, 0.8062122] , [0.0063055, 1.3689392] , [0.0089208, 0.4310608] , [0.0054884, 1.7440906] ,
                    # [0.0102489, 0.0559094] , [0.0048608, 2.0723481] , [0.0095719, 1.1299178] , [0.0104473, 0.9500822] , [0.0080349, 1.4895889] ,
                    # [0.0124457, 0.5904111] , [0.0067447, 1.8492601] , [0.0148264, 0.2307399] , [0.0057870, 2.1639723] , [0.0141035, 1.3563747] ,
                    # [0.0155990, 1.1836253] , [0.0209188, 1.5863747] , [0.0231370, 1.4136253] , [0.0115288, 1.7018736] , [0.0190826, 0.8381264] ,
                    # [0.0171000, 1.9318736] , [0.0283041, 1.0681264] , [0.0094242, 2.0473725] , [0.0233442, 0.4926275] , [0.0139783, 2.2773725] ,
                    # [0.0346250, 0.7226275] , [0.0079003, 2.3496840] , [0.0278469, 0.1903160] , [0.0117181, 2.5796840] , [0.0413036, 0.4203160] ,
                    # [0.0061407, 2.7815576] , [0.0091082, 3.0115576] , [0.0047730, 3.2134311] , [0.0070795, 3.4434311] , [0.0289718, 1.8301501] ,
                    # [0.0326525, 1.6698499] , [0.0405040, 2.0801501] , [0.0456498, 1.9198499] , [0.0228083, 2.1507503] , [0.0414761, 1.3492497] ,
                    # [0.0318871, 2.4007503] , [0.0579857, 1.5992497] , [0.0179560, 2.4713505] , [0.0526842, 1.0286495] , [0.0736552, 1.2786495] ,
                    # [0.0145651, 2.7518757] , [0.0649497, 0.7481243] , [0.0908030, 0.9981243] , [0.0108009, 3.1526260] , [0.0875851, 0.3473740] ,
                    # [0.1224484, 0.5973740] , [0.0491940, 2.3144969] , [0.0573179, 2.1855031] , [0.0607498, 2.5644969] , [0.0707821, 2.4355031] ,
                    # [0.0750202, 2.8144969] , [0.0874090, 2.6855031] , [0.0926426, 3.0644969] , [0.1079417, 2.9355031] , [0.0362373, 2.5724846] ,
                    # [0.0778121, 1.9275154] , [0.0447495, 2.8224846] , [0.0960904, 2.1775154] , [0.0552613, 3.0724846] , [0.1186623, 2.4275154] ,
                    # [0.0682424, 3.3224846] , [0.1056340, 1.6695277] , [0.0329634, 3.0804723] , [0.1304478, 1.9195277] , [0.0407066, 3.3304723] ,
                    # [0.1380276, 1.4437885] , [0.0252272, 3.3062115] , [0.0043682, 2.5774243] , [0.0035859, 3.0366883] , [0.0258983, 2.7759114] ,
                    # [0.0203956, 3.0290436] , [0.0145206, 3.3906610] , [0.0017250, 0.0500000] , [0.0500000, 0.1000000] , [0.0010000, 0.3450000] ,
                    # [0.0025000, 0.5500000] , [0.0050000, 0.7500000] , [0.0075000, 0.9000000] , [0.0100000, 1.0400000] , [0.0220000, 1.5000000] ,
                    # [0.0430000, 2.0000000] , [0.1000000, 3.0000000] ]

# SIM_CONSTS['MUT']['contact inhibition']  = mut_k_lam_vals[{{ k_lam_idx }}][0]
# SIM_CONSTS['MUT']['Cell Stiffness Mean'] = mut_k_lam_vals[{{ k_lam_idx }}][1]



# soft_supercomp_k_lams =[[0.0010221, 0.2474125], [0.0009784, 0.4425875], [0.001547 , 0.5450875], [0.001616 , 0.3499125], 
                        # [0.002446 , 0.6475875], [0.0025552, 0.4524125], [0.0008966, 0.8329374], [0.0014176, 0.9354374], 
                        # [0.0022415, 1.0379374], [0.0027884, 0.0620626], [0.0034389, 0.7460804], [0.0036349, 0.5539196], 
                        # [0.0048633, 0.8460804], [0.0051406, 0.6539196], [0.0030779, 1.1304019], [0.0040612, 0.1695981], 
                        # [0.0057434, 0.2695981], [0.0072442, 0.9937878], [0.0077648, 0.8062122], [0.0089208, 0.4310608], 
                        # [0.0102489, 0.0559094], [0.0095719, 1.1299178], [0.0104473, 0.9500822], [0.001725 , 0.05     ], 
                        # [0.001    , 0.345    ], [0.0025   , 0.55     ], [0.005    , 0.75     ], [0.0075   , 0.9      ], 
                        # [0.01     , 1.04     ], [0.0055   , 1.175    ], [0.00085  , 1.1      ], [0.0012   , 0.05     ],

# SIM_CONSTS['MUT']['contact inhibition']  = soft_supercomp_k_lams[{{ k_lam_idx }}][0]
# SIM_CONSTS['MUT']['Cell Stiffness Mean'] = soft_supercomp_k_lams[{{ k_lam_idx }}][1]




# {{Uncrowded_Growth_Rate_Multiplier}}  {{Sizer_Thresh_Multiplier}}  set in ParameterScanSpecs.json
SIM_CONSTS['MUT']['growth rate']     = SIM_CONSTS['WT']['growth rate']     * {{ Uncrowded_Growth_Rate_Multiplier }}
SIM_CONSTS['MUT']['sizer threshold'] = SIM_CONSTS['WT']['sizer threshold'] * {{          Sizer_Thresh_Multiplier }}
SIM_CONSTS['MUT']['sizer threshold noise'] = 0.025 * SIM_CONSTS['MUT']['sizer threshold']

# Mutant colony size {{init_radius}} set in ParameterScanSpecs.json
SIM_CONSTS['SIM']['START SIM THRESHOLD RADIUS'] = {{ init_radius }}


######################################################################
######            Globals for some plotting classes           ########
######################################################################
global TOTAL_DIV_AND_APO_COUNTS_editable
TOTAL_DIV_AND_APO_COUNTS_editable = { 'WT': { 'DIV': 0, 'APO': 0 } ,
                                     'MUT': { 'DIV': 0, 'APO': 0 } }


global DENSITY_DEATH_RUNNING
DENSITY_DEATH_RUNNING = False

######################################################################
##########           Func Defs for saving data        ################
######################################################################
def make_new_csv_file(filename, col_names):
    try:
        with open(filename, 'w') as file_obj:
            writer_obj = csv.writer(file_obj)
            writer_obj.writerow( col_names )
            return False
    except IOError:
        print("I/O error in make new csv "+filename)
        return True

def update_csv_file(filename, data):
    try:
        with open(filename, 'a+') as file_obj:
            writer_obj = csv.writer(file_obj)
            for line in data:
                writer_obj.writerow( line )
    except IOError:
        print("Failed to update csv "+filename)

def kill_simulation(self, reason_for_murder):
    for i in range(10):
        print(reason_for_murder)
        time.sleep(1)
    self.stopSimulation()  
    

def dict_writer(data, file_path, func_here_flag=False):
    if func_here_flag: print("Entering dict_writer()")

    original_stdout = sys.stdout # Save a reference to the original standard output
    with open(file_path, 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print(data)
        sys.stdout = original_stdout # Reset the standard output to its original value

def dict_reader(file_path):
    with open(file_path, 'r') as f:
        s = f.read()
        file_dic = ast.literal_eval(s)
    return file_dic


def get_subdirectories(dir_path, excluded_folders=[], help=False):
    dir_items = os.listdir( dir_path ) 
    subdirs = []
    for its in dir_items:
        if os.path.isdir( os.path.join(dir_path, its ) ):
            if its not in excluded_folders:
                subdirs.append( its + '/' )

    return subdirs



def get_sim_param_differences(dic1):
    differences = {}
    for k in dic1['WT']:
        if dic1['WT'][k] != dic1['MUT'][k]:
            no_differences = False
            k_str = k.replace('Cell Stiffness Mean','lambda').replace('elastic modulous','lambda')
            k_str = k_str.replace('contact inhibition','k')
            k_str = k_str.replace('sizer threshold','As').replace('growth rate','G')
            if k == 'growth rate':
                differences[k_str] = dic1['MUT'][k]/dic1['WT'][k]
            else:
                differences[k_str] = dic1['MUT'][k]
                
    return differences
    
    
def get_value_from_dist(mean=1, std=0, skew=0, lower=0, upper=2 ):    
    counter = 0
    draw_val_from_dist = True
    while draw_val_from_dist:
        if counter >= 100:
            return 'Error in get_value_from_skew_dist()'
        counter += 1
        
        #draw = spy.stats.skewnorm.rvs(skew, loc=mean, scale=std)
        draw = random.normalvariate( mean, std )
        if (lower <= draw) and (draw <= upper):
            return draw
                    



## Record Values for this simulation before it starts
if SIM_CONSTS['SIM']['SAVE_DATA']:
    dict_writer(SIM_CONSTS, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'Simulation_Parameter_Values_init.txt')   

 
######################################################################
######################################################################
class Initializer(SteppableBasePy):
    def __init__(   self,frequency=1, 
                    piff_start=False, 
                    init_type=1):
        
        SteppableBasePy.__init__(self,frequency)   
        self.piff_start  = piff_start
        self.init_type   = init_type


    def start_sim_from_piff_radius_determines_cell_type(self):
    
        for cell in self.cell_list:
            X, Y = abs(0.5*SIM_CONSTS['SIM']['X_MAX'] - cell.xCOM), abs(0.5*SIM_CONSTS['SIM']['Y_MAX'] - cell.yCOM)
            R = np.sqrt( X**2 + Y**2 )

            CELL_ATTR_TRACKER[str(cell.id)] = {}

            if R <= SIM_CONSTS['SIM']['START SIM THRESHOLD RADIUS']:
                c_species = 'MUT' 
                cell.type = 5
            else:
                c_species = 'WT' 
                cell.type = 1
                
                
            CELL_ATTR_TRACKER[str(cell.id)]['SPECIES'] = c_species            
            CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] = True
            cell.targetVolume = cell.volume 

            skew  = SIM_CONSTS[c_species]['Cell Stiffness Skew']
            mean  = SIM_CONSTS[c_species]['Cell Stiffness Mean']
            std   = SIM_CONSTS[c_species]['Cell Stiffness std' ]
            lower = SIM_CONSTS[c_species]['Cell Stiffness Lower Bound']
            upper = SIM_CONSTS[c_species]['Cell Stiffness Upper Bound']
            #cell.lambdaVolume = get_value_from_dist(mean=mean, std=std, skew=skew, lower=lower, upper=upper )
            cell.lambdaVolume = mean
            #cell.lambdaVolume = SIM_CONSTS[c_species]['elastic modulous']  
            

    
    def start(self): 
    
        self.start_sim_from_piff_radius_determines_cell_type()
    

    def finish(self):
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            dict_writer(SIM_CONSTS, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'Simulation_Parameter_Values_end.txt')    
            


######################################################################
######################################################################
class Simulation_Killer(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    
    def step(self, mcs):

        if KILL_SIM_IF_MUT_ELIMINATED or KILL_SIM_IF_MUT_EXCEEDS_WT:
            cell_counts = {
                         'WT_count': 0,
                        'MUT_count': 0,
                        }

            for cell in self.cell_list: 
                c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES'] 

                cell_counts[c_species + '_count'] += 1

            if KILL_SIM_IF_MUT_ELIMINATED and (cell_counts['MUT_count'] == 0):
                death_message = {'reason for killing': 'mutants eliminated'}
                if SIM_CONSTS['SIM']['SAVE_DATA']:
                    dict_writer(death_message, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'reason_sim_was_killed.txt') 
                kill_simulation(self, 'Mutants eliminated')

            if KILL_SIM_IF_WT_ELIMINATED and (cell_counts['WT_count'] == 0):
                death_message = {'reason for killing': 'wild type eliminated'}
                if SIM_CONSTS['SIM']['SAVE_DATA']:
                    dict_writer(death_message, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'reason_sim_was_killed.txt') 
                kill_simulation(self, 'Wild type eliminated')
                
            if KILL_SIM_IF_MUT_EXCEEDS_WT and (cell_counts['MUT_count'] > cell_counts['WT_count']):
                death_message = {'reason for killing': 'mutants out number wild type cells'}
                if SIM_CONSTS['SIM']['SAVE_DATA']:
                    dict_writer(death_message, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'reason_sim_was_killed.txt') 
                kill_simulation(self, 'Mutants exceed Wild Type')


        if KILL_SIM_IF_R_MUT_EXCEEDS_THRESHOLD:
            mut_R_beyond_thrshold = False

            for cell in self.cell_list: 
                cell_X, cell_Y = abs(0.5*SIM_CONSTS['SIM']['X_MAX'] - cell.xCOM), abs(0.5*SIM_CONSTS['SIM']['Y_MAX'] - cell.yCOM)
                cell_R = np.sqrt( cell_X**2 + cell_Y**2 )

                if cell_R >= R_MUT_THRESHOLD:
                    death_message = {'reason for killing': 'mutants have gone beyond radius threshold'}
                    if SIM_CONSTS['SIM']['SAVE_DATA']:
                        dict_writer(death_message, SIM_CONSTS['SIM']['DATA_OUT_PATH'] + 'reason_sim_was_killed.txt') 
                    kill_simulation(self, 'mutants have gone beyond radius threshold')

        
        
######################################################################
######################################################################
class Crowded_Growth_Steppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        global Ave_Vol_Data
        global Ave_Vol_Data_Filename
        Ave_Vol_Data_Filename = 'Ave_Vol_Data.csv'        
        Ave_Vol_Col_Names = ['Time',  'WT Count',  'WT Ave Vol',  'WT Total Vol',  'WT Ave dA/dt',  'WT Total dA/dt',
                                     'MUT Count', 'MUT Ave Vol', 'MUT Total Vol', 'MUT Ave dA/dt', 'MUT Total dA/dt']
                                     
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Ave_Vol_Data_Filename, Ave_Vol_Col_Names)           
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Ave_Vol_Data_Filename + 'Killing Sim')
        Ave_Vol_Data = []
    
        if SIM_CONSTS['SIM']['SHOW_AVE_AREA_PLOT']:
            self.pW_ave_area = self.add_new_plot_window( title     = 'Average Area',
                                                     x_axis_title = 'Time (mcs)'  ,
                                                     y_axis_title = 'Ave Area'       ,
                                                     x_scale_type = 'linear'        ,
                                                     y_scale_type = 'linear'        )
            
            self.pW_ave_area.add_plot('Avg WT Area' ,  style='Dots',  color='green',  size=5 ) 
            self.pW_ave_area.add_plot('Avg MUT Area',  style='Dots',  color='purple'  ,  size=5 )  

        if SIM_CONSTS['SIM']['SHOW_AVE_ADDED_dATdt_PLOT']:
            self.pW_added_dATdt = self.add_new_plot_window( title   = 'Average Added dAT/dt',
                                                     x_axis_title = 'Time (mcs)',
                                                     y_axis_title = 'Ave Added Area'  ,
                                                     x_scale_type = 'linear'      ,
                                                     y_scale_type = 'linear'      )
            
            self.pW_added_dATdt.add_plot('Avg WT Added dAT/dt' ,  style='Dots',  color='green',  size=5 ) 
            self.pW_added_dATdt.add_plot('Avg MUT Added dAT/dt',  style='Dots',  color='Purple'  ,  size=5 ) 
        
    def step(self,mcs):    
        global Ave_Vol_Data
    
        plot_data = {}
        for cs in ['WT', 'MUT']:
            plot_data[cs] ={'count'           : 0 ,
                            'total_area'      : 0 ,
                            'added_dATdt'     : 0 ,
                            'ave_area'        : 0 ,
                            'ave_added_dATdt' : 0 }
                            
        for cell in self.cell_list:   
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
                c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES'] 
                
                G = random.normalvariate( SIM_CONSTS[c_species]['growth rate'], 0.05*SIM_CONSTS[c_species]['growth rate'] )
                if G < 0:
                    G = 0.
                k = SIM_CONSTS[c_species]['contact inhibition']
                dATdt = G * np.exp( -k * (cell.volume - cell.targetVolume)**2 )
                cell.targetVolume += dATdt
                
                plot_data[c_species]['added_dATdt'] += dATdt
                plot_data[c_species]['count']       += 1
                plot_data[c_species]['total_area']  += cell.volume


        ## END SIMULATION IF ALL CELLS ARE DEAD
        if (plot_data['WT']['count'] + plot_data['MUT']['count']) == 0 and SIM_CONSTS['SIM']['END_SIM_IF_NO_CELLS']:
            kill_simulation(self, 'There are no cells left in the simulation...  Killing Sim')


        for cs in ['WT', 'MUT']:
            if plot_data[cs]['count'] != 0:
                plot_data[cs]['ave_area']        = plot_data[cs]['total_area']  / float( plot_data[cs]['count'] ) 
                plot_data[cs]['ave_added_dATdt'] = plot_data[cs]['added_dATdt'] / float( plot_data[cs]['count'] ) 
                if SIM_CONSTS['SIM']['SHOW_AVE_AREA_PLOT']:     
                    self.pW_ave_area.add_data_point(  'Avg '+cs+' Area',  mcs, plot_data[cs]['ave_area'])  
                if SIM_CONSTS['SIM']['SHOW_AVE_ADDED_dATdt_PLOT']:     
                    self.pW_added_dATdt.add_data_point(  'Avg '+cs+' Added dAT/dt',  mcs, plot_data[cs]['ave_added_dATdt'])  
        
                     
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            Ave_Vol_Data.append( [ mcs,  plot_data['WT']['count'],  plot_data['WT']['ave_area'],  plot_data['WT']['total_area'],  plot_data['WT']['ave_added_dATdt'],  plot_data['WT']['added_dATdt'],   
                                        plot_data['MUT']['count'], plot_data['MUT']['ave_area'], plot_data['MUT']['total_area'], plot_data['MUT']['ave_added_dATdt'], plot_data['MUT']['added_dATdt'],    ] )        
            
            if mcs % SIM_CONSTS['SIM']['SAVE_DATA_TO_CSV_INTERVAL']:
                update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Ave_Vol_Data_Filename, Ave_Vol_Data )
                Ave_Vol_Data = []
            
        


            
######################################################################
######################################################################
class Sizer_Timer_Mitosis_Steppable(MitosisSteppableBase):  
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)
        
    def get_distance_to_center(self, x, y):
        return np.sqrt( (SIM_CONSTS['SIM']['X_MAX']*.5 - x)**2 + (SIM_CONSTS['SIM']['Y_MAX']*.5 - y)**2 )
        
    def get_distance_between_points(self, x1, y1, x2, y2):
        return np.sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
        
    def start(self):   
        ##Initialize data save files 
        global Cell_Heritage_Data 
        global Heritage_Data_Filename
        Heritage_Data_Filename = 'Cell_Heritage_Data__Timer_Sizer_Mitosis.csv'
        Cell_Heritage_Col_Names = ['Time', 
                                        'Mother ID', 'Mother Species', 'Mother Type', 'Mother Volume', 'Mother Target Volume', 
                                                     'Mother X'      , 'Mother Y'   , 'Mother R'     , 
                                                     'Mother Sizer Threshold'       , 'Mother Timer Threshold', 'Mother Start of Timer Phase' ,
                                                     'Mother Stiffness',
                                                     
                                        'Parent ID', 'Parent Species', 'Parent Type', 'Parent Volume', 'Parent Target Volume', 
                                                     'Parent X',       'Parent Y'   , 'Parent R', 
                                                     'Parent Sizer Threshold'       , 'Parent Timer Threshold',  
                                                     'Parent Stiffness',
                                                     
                                        'Child ID' , 'Child Species', 'Child Type',  'Child Volume',  'Child Target Volume',  
                                                     'Child X',     'Child Y',       'Child R', 
                                                     'Child Sizer Threshold',        'Child Timer Threshold',
                                                     'Child Stiffness'
                                    ]
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Heritage_Data_Filename, Cell_Heritage_Col_Names)      
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Heritage_Data_Filename + 'Killing Sim')
        Cell_Heritage_Data = []            
        
        global Cell_Cycle_Data_Filename 
        global Cell_Cycle_Data
        global Cell_Cycle_Col_Names
        Cell_Cycle_Data_Filename = 'Cell_Cycle_Data.csv'
        Cell_Cycle_Col_Names = [ 'Cycle Time'  , 'Cycle Area'   , 'Cycle Movement',              'Cycle R', 
                                 'Birth mcs'   , 'Birth Area'   , 'Birth X'   , 'Birth Y'   ,    'Birth R',
                                 'Division mcs', 'Division Area', 'Division X', 'Division Y', 'Division R',
                                 'Sizer Threshold', 
                                 'Timer Threshold',                                   
                                 'Start of Timer Phase',
                                 'Cell Stiffness' ] 
                                 
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Cell_Cycle_Data_Filename, Cell_Cycle_Col_Names)      
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Cell_Cycle_Data_Filename + 'Killing Sim') 
        Cell_Cycle_Data = []
        
        
        
        ## Set sizer-timer thresholds for each cell at start of sim
        ## Set initial cell's values in CELL_CYCLE_DATA_DICT for data saving
        for cell in self.cell_list:
            c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES'] 
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(cell.id)] = { col:None for col in Cell_Cycle_Col_Names }
            
            
            ### This is way overly complicated... sorry future me... 
            has_sizer_threshold = str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds']
            has_timer_threshold = str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds']
            
            if has_sizer_threshold:
                sizer_threshold = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( cell.id )] 
            else:
                sizer_threshold = int( random.normalvariate(SIM_CONSTS[c_species]['sizer threshold'], SIM_CONSTS[c_species]['sizer threshold noise']) )
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( cell.id )]   = sizer_threshold
                
            
            if has_timer_threshold:
                timer_threshold = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( cell.id )] 
            else:
                timer_threshold = int( random.normalvariate(SIM_CONSTS[c_species]['timer threshold'], SIM_CONSTS[c_species]['timer threshold noise']) )
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( cell.id )]   = timer_threshold
            
            
            
        
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( cell.id )]['Sizer Threshold'] = sizer_threshold
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( cell.id )]['Timer Threshold'] = timer_threshold

            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( cell.id )]['Cell Stiffness'] = cell.lambdaVolume
            
            
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth mcs']  = 0  #### OTHER HARDCODED COSNT
            ################ Target Vol because the cell is a small square at mcs=0
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth Area'] = cell.targetVolume  
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth X']    = cell.xCOM
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth Y']    = cell.yCOM
            CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( cell.id ) ]['Birth R']    = self.get_distance_to_center( cell.xCOM, cell.yCOM )
                                 
        
        ## Initialize plots
        global total_mitosis_count
        total_mitosis_count = { 'WT':0, 'MUT':0 }      
        
        if SIM_CONSTS['SIM']['SHOW_DIVISION_COUNT_PLOT']:
            self.pW_div_count = self.add_new_plot_window(   title      = 'Cell Division Count',
                                                         x_axis_title = 'Time (mcs)'         ,
                                                         y_axis_title = 'Division Count'     ,
                                                         x_scale_type = 'linear'             ,
                                                         y_scale_type = 'linear'             )            
            self.pW_div_count.add_plot('WT Mitosis'    ,  style='Dots',  color='green',  size=5 ) 
            self.pW_div_count.add_plot('MUT Mitosis',  style='Dots',  color='purple'  ,  size=5 ) 

        if SIM_CONSTS['SIM']['SHOW_POPULATION_COUNT_PLOT']:
            self.pW_pop_count = self.add_new_plot_window( title   = 'Population Count',
                                                     x_axis_title = 'Time (mcs)'      ,
                                                     y_axis_title = 'Population Count',
                                                     x_scale_type = 'linear'          ,
                                                     y_scale_type = 'linear'          )            
            self.pW_pop_count.add_plot('WT Count'    ,  style='Dots',  color='green',  size=5 ) 
            self.pW_pop_count.add_plot('MUT Count',  style='Dots',  color='purple'  ,  size=5 ) 
        
        
    def step(self,mcs):
        global total_mitosis_count
        global TOTAL_DIV_AND_APO_COUNTS_editable
        global Cell_Heritage_Data 
        global Cell_Cycle_Data         
        global Cell_Cycle_Col_Names
        
        global CELL_ATTR_TRACKER
        
        plot_data = {}
        for cs in ['WT', 'MUT']:
            plot_data[cs] = { 'div_count': 0 ,
                                  'count': 0 }
        cells_to_divide = []  

        ## Check if Mitosis Conditions Met for each cell
        for cell in self.cell_list:
            c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES'] 
            plot_data[c_species]['count'] += 1
            
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
                
                ## Logic for Sizer-Timer Mitosis
                if str( cell.id ) in  CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs']:
                    time_past = mcs - CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]
                    time_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str(cell.id)]
                    if time_past >= time_thresh:
                        cells_to_divide.append( cell )
                        CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'].pop( str(cell.id), None )
                        
                        plot_data[c_species]['div_count'] += 1
                else:
                    if cell.volume > CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]:
                        # print("WT bigger than min area - mcs:",mcs)
                        CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]         = mcs
                        CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(cell.id)]['Start of Timer Phase'] = mcs
                        # print("dict:",MITOSIS_TIMER_DICT)
                    
        
        
        for cs in ['WT', 'MUT']:
            TOTAL_DIV_AND_APO_COUNTS_editable[cs]['DIV'] +=  plot_data[cs]['div_count']
            
        
        for cell in cells_to_divide:       
        
            mother_id         = cell.id
            mother_species    = CELL_ATTR_TRACKER[str(mother_id)]['SPECIES'] 
            mother_type       = cell.type         
            mother_stiffness  = cell.lambdaVolume
            mother_vol        = cell.volume
            mother_target_vol = cell.targetVolume
            mother_x          = cell.xCOM
            mother_y          = cell.yCOM
            mother_R          = self.get_distance_to_center( mother_x, mother_y )
            mother_sizer_threshold   = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Sizer Threshold']
            mother_timer_threshold   = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Timer Threshold']
            mother_start_timer_phase = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( mother_id )]['Start of Timer Phase']
        
            self.divideCellAlongMinorAxis(cell) #No calling the "cell" object after this            
            
            ############################################################################
            ## NOTE: "parentCell" is a terrible name because it is a daughter cell.    
            ##       CC3D gives the mother cell id to the "parentCell" and generates 
            ##       a new unique cell id for the "childCell".
            ############################################################################
            self.parentCell.targetVolume /= 2.     #### OTHER HARDCODED COSNT
            self.cloneParent2Child()  
            
            self.childCell.type = self.parentCell.type
            parent_id = self.parentCell.id
            child_id  = self.childCell.id
            birth_target_vol = self.childCell.targetVolume   
            
            #### OTHER HARDCODED COSNT
            if CELL_ATTR_TRACKER['COLOR_CELLS'] == {}:
                if self.parentCell.type in [4,8]:
                    self.parentCell.type -= 3
                else:
                    self.parentCell.type += 1
            else:
                #CELL_ATTR_TRACKER['COLOR_CELLS'][str( self.childCell.id )] = [ x for x in CELL_ATTR_TRACKER['COLOR_CELLS'][str( parent_id )] ]       
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(  self.childCell.id )] = [ self.childCell.volume  for x in range(9) ]         
                CELL_ATTR_TRACKER['COLOR_CELLS'][str( self.parentCell.id )] = [ self.parentCell.volume for x in range(9) ]       
                    
             
            
            # only need to set child because   mother_id = parent_id
            CELL_ATTR_TRACKER[str(child_id)] = {} 
            CELL_ATTR_TRACKER[str(child_id)]['SPECIES'] = mother_species                     
            CELL_ATTR_TRACKER[str(child_id)]['IS_ALIVE'] = True
            
            
            ################################################################################
            # 'Cycle Time'  ,   'Cycle  Area',      'Cycle Movement'     ,    'Cycle R',   #
            #                                                                              #
            # 'Birth mcs'   ,    'Birth Area',    'Birth X',    'Birth Y',    'Birth R',   #
            #                                                                              #
            # 'Division mcs', 'Division Area', 'Division X', 'Division Y', 'Division R',   #
            #                                                                              #
            # 'Sizer Threshold',                                                           #
            # 'Timer Threshold',                                                           #
            # 'Start of Timer Phase'                                                       #
            ################################################################################
            if str(cell.id) in CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT']:            
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division mcs']  = mcs
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division Area'] = mother_vol
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cell Stiffness']= mother_stiffness
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division X']    = mother_x
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division Y']    = mother_y
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division R']    = mother_R
                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Time']     = mcs - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth mcs']
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Area']     = mother_vol - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth Area']

                
                
                birth_x, birth_y = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth X'] , CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth Y']                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle Movement'] = self.get_distance_between_points( mother_x, mother_y, birth_x, birth_y ) 
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Cycle R']        = CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Division R'] - CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)]['Birth R'] 
                
                ## Store cell's cycle data in a temp list
                c_cyc_data = [ CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)][col] for col in Cell_Cycle_Col_Names ]
                ## Check that there are no "None" values in temp list
                if None in c_cyc_data:
                    None_cols = [ col for col in Cell_Cycle_Col_Names if CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(mother_id)][col] == None ]
                    err_msg = 'None in Cell Cycle Data entry columns: '
                    for k in None_cols:
                        err_msg += k + ', '
                    kill_simulation(self, err_msg) # End sim if there are "None" values
                else:
                    Cell_Cycle_Data.append( c_cyc_data ) # Save the cycle data if all is kosher
                
            ## Initialize Cycle data for daughter cells
            for divided_cell in [ self.parentCell, self.childCell ]:               
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str(divided_cell.id)] = { col:None for col in Cell_Cycle_Col_Names }
            
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth mcs']  = mcs
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth Area'] = divided_cell.volume
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth X']    = divided_cell.xCOM
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth Y']    = divided_cell.yCOM
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][ str( divided_cell.id ) ]['Birth R']    = self.get_distance_to_center( divided_cell.xCOM, divided_cell.yCOM )
           
                
                
                sizer_threshold = int( random.normalvariate(SIM_CONSTS[mother_species]['sizer threshold'], SIM_CONSTS[mother_species]['sizer threshold noise']) )
                timer_threshold = int( random.normalvariate(SIM_CONSTS[mother_species]['timer threshold'], SIM_CONSTS[mother_species]['timer threshold noise']) )
                    
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( divided_cell.id )]   = sizer_threshold
                CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( divided_cell.id )]   = timer_threshold
                
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( divided_cell.id )]['Sizer Threshold'] = sizer_threshold
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( divided_cell.id )]['Timer Threshold'] = timer_threshold        
            
                
                skew  = SIM_CONSTS[mother_species]['Cell Stiffness Skew']
                mean  = SIM_CONSTS[mother_species]['Cell Stiffness Mean']
                std   = SIM_CONSTS[mother_species]['Cell Stiffness std' ]
                lower = SIM_CONSTS[mother_species]['Cell Stiffness Lower Bound']
                upper = SIM_CONSTS[mother_species]['Cell Stiffness Upper Bound']
                #divided_cell.lambdaVolume = get_value_from_dist(mean=mean, std=std, skew=skew, lower=lower, upper=upper )
                divided_cell.lambdaVolume = mean
                CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( divided_cell.id )]['Cell Stiffness'] = divided_cell.lambdaVolume

            
            Cell_Heritage_Data.append( [mcs, 
                mother_id,  mother_species, mother_type,  mother_vol, mother_target_vol  ,          
                    mother_x   ,  mother_y  , mother_R, 
                    mother_sizer_threshold, mother_timer_threshold, mother_start_timer_phase,
                    mother_stiffness,
                            
                parent_id, mother_species, self.parentCell.type, self.parentCell.volume, self.parentCell.targetVolume, 
                    self.parentCell.xCOM, self.parentCell.yCOM, CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( parent_id )]['Birth R'],
                    CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( parent_id )], CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( parent_id )],
                    CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( parent_id )]['Cell Stiffness'],
                            
                child_id, mother_species, self.childCell.type,  self.childCell.volume,  self.childCell.targetVolume,  
                    self.childCell.xCOM,  self.childCell.yCOM, CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( child_id )]['Birth R'],
                    CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str( child_id )], CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str( child_id )],
                    CELL_ATTR_TRACKER['CELL_CYCLE_DATA_DICT'][str( child_id )]['Cell Stiffness'],
                                        ] )
        
            
        ## Live plots
        for cs in ['WT', 'MUT']:
            if SIM_CONSTS['SIM']['SHOW_DIVISION_COUNT_PLOT'] and (plot_data[cs]['count'] != 0):
                total_mitosis_count[cs]  += plot_data[cs]['div_count']       
                self.pW_div_count.add_data_point( cs+' Mitosis'    , mcs, total_mitosis_count[cs] ) 
                
                
                
            if SIM_CONSTS['SIM']['SHOW_POPULATION_COUNT_PLOT'] and (plot_data[cs]['count'] != 0):
                self.pW_pop_count.add_data_point( cs+' Count'    , mcs, plot_data[cs]['count']  )
        
        ## Save data
        if SIM_CONSTS['SIM']['SAVE_DATA'] and (mcs % SIM_CONSTS['SIM']['SAVE_DATA_TO_CSV_INTERVAL']):
            update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Heritage_Data_Filename, Cell_Heritage_Data )            
            Cell_Heritage_Data = []
            
            update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Cell_Cycle_Data_Filename, Cell_Cycle_Data )            
            Cell_Cycle_Data = []
                
            


        
###############################################################################
###############################################################################
###############################################################################
class Death_Steppable(SteppableBasePy):
    def __init__(self,  frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.calc_prob_of_contact_death_freq =  frequency
        
        ## fast cycle times for MDCK cells is ~10 hr = 600 min
        min_per_mcs = 600 / SIM_CONSTS['WT']['timer threshold']
        ## This gets the number of times P_mech( rho ) needs to check 
        #       if the cell dies. This is because  1 mcs ~= 2.6 min
        #       and the experimental data we fit is given in terms
        #       of a probability every 4 minutes.  So, with freq = 10 mcs
        #       then each cell will need to check if it dies 6 times
        #       when the Death_Steppable() is called.
        self.check_prob_of_mech_death_X_times = int( frequency * min_per_mcs // SIM_CONSTS['SIM']['P_MECH_MIN_PER_FRAME'] )
        
        
    def start(self):
        #########################################
        global Local_Density_Data
        global Local_Density_Data_Filename
        Local_Density_Data_Filename = 'Local_Density_Data.csv'        
        Local_Density_Col_Names = ['Time', 'ID', 'Type', 'Species',
                     'Local Density', 'Number of Neighbors', 'Neighbor IDs',
                     'WT Neighbors', 'MUT Neighbors', 
                     'WT Contact Length', 'MUT Contact Length', 
                     'Medium Contact Length', 'Perimeter Length',
                     'Volume', 'Target Volume', 'Cell Stiffness', 'Stress', 'Pressure', 
                     'X', 'Y', 'R' ]
                     
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Local_Density_Data_Filename, Local_Density_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Local_Density_Data_Filename + 'Killing Sim')  
        Local_Density_Data = []
        
        
        #######################################
        global P_Apo_Density_Data
        global P_Apo_Density_Data_Filename
        P_Apo_Density_Data = []
        P_Apo_Density_Data_Filename = 'P_Apo_Density_Data.csv'  

        P_Apo_Density_Col_Names = Local_Density_Col_Names            
        if SIM_CONSTS['SIM']['SAVE_DATA'] and SIM_CONSTS['SIM']['CALCULATE_MECH_PROB_APO']:      
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + P_Apo_Density_Data_Filename, P_Apo_Density_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + P_Apo_Density_Data_Filename + 'Killing Sim')  
            
            
        
        global density_death_list  
        global death_time_dict
        density_death_list = []
        death_time_dict    = {}
        
        global tot_wt_apo_den_count_plot 
        global tot_mut_apo_den_count_plot  
        tot_wt_apo_den_count_plot, tot_mut_apo_den_count_plot = 0, 0
        
        global DENSITY_DEATH_RUNNING
        DENSITY_DEATH_RUNNING = True
        
        


        #########################################
        global Contact_Death_Data
        Contact_Death_Data = []
        global Contact_Death_Data_Filename
        Contact_Death_Data_Filename = 'P_Apo_Contact_Data.csv'  

        Contact_Death_Col_Names = [ 'Time', 'ID', 'Species',
                                    'Volume', 'Target Volume', 'Cell Stiffness', 
                                    'Number of Neighbors', 'Neighbor IDs', 'WT Neighbors',    'MUT Neighbors',
                                    'Total Perimeter length', 'Medium Contact Area', 'WT Contact Area', 'MUT Contact Area',
                                    'X', 'Y', 'R' ] 

        if SIM_CONSTS['SIM']['SAVE_DATA'] and SIM_CONSTS['SIM']['CALCULATE_BIOCHEM_PROB_APO']:           
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Contact_Death_Data_Filename, Contact_Death_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Contact_Death_Data_Filename + 'Killing Sim')

        global total_wt_contact_apo_count
        global total_mut_contact_apo_count
        total_wt_contact_apo_count  = 0
        total_mut_contact_apo_count = 0



        
        
        #########################################         
        global Vol_Den_Count_Data
        global Vol_Den_Count_Data_Filename
        Vol_Den_Count_Data_Filename = 'Size_Volume_Count_Data.csv'        
        Col_Names = ['Time', 'Average WT Density', 'Average WT Volume', 'WT Count', 'Average MUT Density', 'Average MUT Volume', 'MUT Count']
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Vol_Den_Count_Data_Filename, Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Vol_Den_Count_Data_Filename + 'Killing Sim')  
        Vol_Den_Count_Data = []
        
        
        
        #########################################
        if SIM_CONSTS['SIM']['SHOW_DENSITY_APO_COUNT_PLOT']:
            self.pW=self.add_new_plot_window( title      = 'Density Apoptosis',
                                           x_axis_title = 'Time (mcs)'       ,
                                           y_axis_title = 'Apoptosis Count'  ,
                                           x_scale_type = 'linear'           ,
                                           y_scale_type = 'linear'           )
            
            self.pW.add_plot('WT Apoptosis'  , style='Dots', color='green', size=5)   
            self.pW.add_plot('MUT Apoptosis', style='Dots', color='purple', size=3) 
            
        if SIM_CONSTS['SIM']['SHOW_AVE_DENSITY_PLOT']:
            self.pW_density=self.add_new_plot_window( title      = 'Average Local Density',
                                                   x_axis_title = 'Time (mcs)'           ,
                                                   y_axis_title = 'Ave Density'          ,
                                                   x_scale_type = 'linear'               ,
                                                   y_scale_type = 'linear'               )       
            self.pW_density.add_plot('WT Density' , style='Dots', color='green' , size=5)   
            self.pW_density.add_plot('MUT Density', style='Dots', color='purple', size=3)   

        if SIM_CONSTS['SIM']['SHOW_CONTACT_APO_COUNT_PLOT']:
            self.pW_biochem=self.add_new_plot_window( title     = 'Contact Apoptosis',
                                           x_axis_title = 'Time (mcs)'      ,
                                           y_axis_title = 'Biochem Apo Count' ,
                                           x_scale_type = 'linear'          ,
                                           y_scale_type = 'linear'           )
            
            
            self.pW_biochem.add_plot( 'WT Biochem Apo',  style='Dots',  color='green' ,  size=5)   
            self.pW_biochem.add_plot('MUT Biochem Apo',  style='Dots',  color='purple',  size=3)       
        
    def death_prob(self, local_den, papo_max, alpha, rho_half):
        return papo_max/(1 + np.exp(-alpha * (local_den - rho_half)))
        
    def get_distance_to_center(self, x, y):
        return np.sqrt( (SIM_CONSTS['SIM']['X_MAX']*.5 - x)**2 + (SIM_CONSTS['SIM']['Y_MAX']*.5 - y)**2 )

    def step(self, mcs):        
        global Local_Density_Data   
        global Local_Density_Data_Filename        
        global Vol_Den_Count_Data
        global Vol_Den_Count_Data_Filename
        global Contact_Death_Data
        global Contact_Death_Data_Filename
        
        global P_Apo_Density_Data
        
        global density_death_list
        global death_time_dict
        
        global tot_wt_apo_den_count_plot 
        global tot_mut_apo_den_count_plot

        global total_wt_contact_apo_count
        global total_mut_contact_apo_count
        
        global TOTAL_DIV_AND_APO_COUNTS_editable
                    
        plot_data = {}
        for cs in ['WT', 'MUT']:
            plot_data[cs] = {     'count': 0 ,
                      'density apo count': 0 ,
                      'contact apo count': 0 ,
                         'summed density': 0 ,
                            'ave density': 0 ,
                              'total vol': 0 ,
                                'ave vol': 0 }
        
        
        # alive_cells_list = [ int(c.id) for c in self.cell_list ]
        # for dying_cell_id in list( death_time_dict.keys() ): 
            # if not( int(dying_cell_id) in alive_cells_list ):
                # death_time_dict[ str( dying_cell_id ) ].insert(0, mcs)
                # P_Apo_Density_Data.append( death_time_dict[ str( dying_cell_id ) ] ) 
                # death_time_dict.pop( str( dying_cell_id ), None)        
        
        for cell in self.cell_list:             
            c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES']        
            
            # This ignores cells that are already dying
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] and cell.volume > 1: #### OTHER HARDCODED COSNT
                
                loc_density = cell.volume**-1.            
                num_neighbors = {'WT neighbor count': 0,
                                'MUT neighbor count': 0,
                              'total neighbor count': 0,
                              
                                   'WT ComSurf Area': 0,
                                  'MUT ComSurf Area': 0,
                               'Medium ComSurf Area': 0,
                               
                               'neighbor IDs':''}
                
                

                for neighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(cell):  
                    ## evaluates True for all cells, but not medium
                    if neighbor:
                        num_neighbors['total neighbor count'] += 1

                        num_neighbors['neighbor IDs'] += str(neighbor.id) + '-'
                        neigh_sp = CELL_ATTR_TRACKER[str(neighbor.id)]['SPECIES']
                        num_neighbors[neigh_sp+' neighbor count'] += 1
                        num_neighbors[neigh_sp+' ComSurf Area'] += commonSurfaceArea

                        if (neighbor.volume > 0) and CELL_ATTR_TRACKER[str(neighbor.id)]['IS_ALIVE']: #### OTHER HARDCODED COSNT
                            loc_density += neighbor.volume**-1.
                    ## Account for the medium medium
                    else:
                        num_neighbors['Medium ComSurf Area'] += commonSurfaceArea

                


                stress = (cell.volume - cell.targetVolume)**2                
                pressure = - 2 * cell.lambdaVolume * (cell.volume - cell.targetVolume)   
                cell_R = self.get_distance_to_center( cell.xCOM, cell.yCOM )
                
                perim_len = sum([num_neighbors[k+' ComSurf Area'] for k in ['WT', 'MUT', 'Medium']])
                
                
                local_den_data_list = [mcs, cell.id, cell.type, c_species,
                                            loc_density, num_neighbors['total neighbor count'], num_neighbors['neighbor IDs'],
                                            num_neighbors['WT neighbor count'], num_neighbors['MUT neighbor count'], 
                                            num_neighbors['WT ComSurf Area'], num_neighbors['MUT ComSurf Area'], 
                                            num_neighbors['Medium ComSurf Area'], perim_len,
                                            cell.volume, cell.targetVolume, cell.lambdaVolume, stress, pressure, 
                                            cell.xCOM, cell.yCOM, cell_R ]
                Local_Density_Data.append( local_den_data_list )
                
                
                       
                plot_data[c_species]['count']          += 1
                plot_data[c_species]['summed density'] += loc_density
                plot_data[c_species]['total vol']      += cell.volume
                
                
                if SIM_CONSTS['SIM']['CALCULATE_MECH_PROB_APO']:
                    #### OTHER HARDCODED COSNTS
                    in_X_bool = (cell.xCOM > SIM_CONSTS['SIM']['EDGE_BUFFER'] and cell.xCOM < SIM_CONSTS['SIM']['X_MAX']-SIM_CONSTS['SIM']['EDGE_BUFFER'])
                    in_Y_bool = (cell.yCOM > SIM_CONSTS['SIM']['EDGE_BUFFER'] and cell.yCOM < SIM_CONSTS['SIM']['Y_MAX']-SIM_CONSTS['SIM']['EDGE_BUFFER'])  
                    if in_X_bool and in_Y_bool:     
                        pmax     = SIM_CONSTS[c_species]['density p_apo,max']
                        alpha    = SIM_CONSTS[c_species]['density death sensitivity']
                        rho_Half = SIM_CONSTS[c_species]['density death susceptibility']
                        
                        #################################################################################
                        ### Need this to match how the probability of death was measured in experiment  #
                        #################################################################################
                        cell_killed = False                                                             #
                        for chck in range( self.check_prob_of_mech_death_X_times ):                     #  
                            if self.death_prob(loc_density, pmax, alpha, rho_Half) > random.random():   #
                                cell_killed = True                                                      #                            
                        if cell_killed:                                                                 # 
                            CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] = False                         # 
                            cell.targetVolume = 0                                                       #
                            cell.lambdaVolume = 5                                                       #     
                            plot_data[c_species]['density apo count'] += 1                              #
                        #################################################################################
                            
                            if CELL_ATTR_TRACKER['COLOR_CELLS'] != {}:
                                cell.type = 9
                                
                        if not( CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] ):                                              
                                density_death_list.append( [ mcs, cell.id ] )                             
                                P_Apo_Density_Data.append( local_den_data_list )
                            
                            

                            
                if SIM_CONSTS['SIM']['CALCULATE_BIOCHEM_PROB_APO']:
                    perimeter_total = sum([ num_neighbors[k] for k in ['Medium ComSurf Area', 'MUT ComSurf Area', 'WT ComSurf Area' ] ]) *1.
                    safe_contact = num_neighbors['Medium ComSurf Area']   +   num_neighbors[c_species+' ComSurf Area'] *1.
                    
                    perimeter_percentage = 1.0 - ( safe_contact / perimeter_total )
                    
                    pmx   = SIM_CONSTS[c_species]['biochemical p_apo,max']
                    nH    = SIM_CONSTS[c_species]['biochemical Hill coef']
                    stp_n = SIM_CONSTS[c_species]['biochemical steepness']**nH
                    p_n = perimeter_percentage**nH
                    Hill_F = pmx * p_n /(stp_n + p_n)
                    
                    time_thresh = SIM_CONSTS[c_species]['timer threshold'] *1.
                    probability_time_fraction = self.calc_prob_of_contact_death_freq / time_thresh
                    if ( probability_time_fraction * Hill_F ) > random.random():
                        CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] = False
                        cell.targetVolume = 0
                        cell.lambdaVolume = 5
                        plot_data[c_species]['contact apo count'] += 1

                        
                        if CELL_ATTR_TRACKER['COLOR_CELLS'] != {}:
                            cell.type = 10

                                        
                        tot_nn = num_neighbors['WT neighbor count'] + num_neighbors['MUT neighbor count']
                        Contact_Death_Data.append([ mcs, cell.id, c_species,
                                                    cell.volume, cell.targetVolume, cell.lambdaVolume,
                                                    tot_nn, num_neighbors['neighbor IDs'],
                                                    num_neighbors['WT neighbor count'], num_neighbors['MUT neighbor count'],
                                                    perimeter_total, num_neighbors['Medium ComSurf Area'], num_neighbors['WT ComSurf Area'], num_neighbors['MUT ComSurf Area'],
                                                    cell.xCOM, cell.yCOM, cell_R  ])
                            
        for cs in ['WT', 'MUT']:
            TOTAL_DIV_AND_APO_COUNTS_editable[cs]['APO'] +=  plot_data[cs]['density apo count']
            TOTAL_DIV_AND_APO_COUNTS_editable[cs]['APO'] +=  plot_data[cs]['contact apo count']
            
            if plot_data[cs]['count'] != 0:
                plot_data[cs]['ave density'] = plot_data[cs]['summed density'] / float( plot_data[cs]['count'] ) 
                plot_data[cs]['ave vol'] = plot_data[cs]['total vol'] / float(plot_data[cs]['count'])
                
                if SIM_CONSTS['SIM']['SHOW_DENSITY_APO_COUNT_PLOT']:
                    self.pW.add_data_point( cs+' Apoptosis'  , mcs, plot_data[cs]['density apo count'] )
                
                if SIM_CONSTS['SIM']['SHOW_AVE_DENSITY_PLOT']:
                    self.pW_density.add_data_point(  cs+' Density' , mcs, plot_data[cs]['ave density']  )  
                
                if SIM_CONSTS['SIM']['SHOW_CONTACT_APO_COUNT_PLOT'] and (plot_data[cs]['count'] != 0):
                    self.pW_biochem.add_data_point(  cs+' Biochem Apo'  , mcs,  plot_data[cs]['contact apo count'])  

        
                            
        Vol_Den_Count_Data.append([ mcs,  plot_data['WT']['ave density'],  plot_data['WT']['ave vol'],  plot_data['WT']['count'], 
                                         plot_data['MUT']['ave density'], plot_data['MUT']['ave vol'], plot_data['MUT']['count']])
        
        if SIM_CONSTS['SIM']['SAVE_DATA'] and (mcs % SIM_CONSTS['SIM']['SAVE_DATA_TO_CSV_INTERVAL']):
            update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Local_Density_Data_Filename, Local_Density_Data )     
            update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Vol_Den_Count_Data_Filename, Vol_Den_Count_Data )    
            
            Local_Density_Data = []
            Vol_Den_Count_Data = []   

            if SIM_CONSTS['SIM']['CALCULATE_MECH_PROB_APO'] and P_Apo_Density_Data != []:
                update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + P_Apo_Density_Data_Filename, P_Apo_Density_Data )
                P_Apo_Density_Data = []      

            if SIM_CONSTS['SIM']['CALCULATE_BIOCHEM_PROB_APO'] and Contact_Death_Data != []:
                update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Contact_Death_Data_Filename, Contact_Death_Data )  
                Contact_Death_Data = []
                 
        
                          

         
            
            
          
       




  


class Extrusion(SteppableBasePy):        
    def __init__(self, frequency=20):
        SteppableBasePy.__init__(self, frequency)
        
        
    def start(self):           
        global Extrusion_Data
        global Extrusion_Data_Filename
        Extrusion_Data_Filename = 'Extrusion_Data.csv'        
        Extrusion_Col_Names = ['Time', 'ID', 'Cell Species', 'Cell Type', 'Volume', 'Average Volume', 'Stress']
        
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Extrusion_Data_Filename, Extrusion_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Extrusion_Data_Filename + 'Killing Sim') 
        Extrusion_Data = []
        
        
    def get_wt_and_mut_average_area(self):
    
        
        plot_data = {}
        for cs in ['WT', 'MUT']:
            plot_data[cs] = {     'count': 0. ,
                               'ave area': 0. ,
                             'total area': 0. }
            
        for cell in self.cell_list:
            c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES']
            plot_data[c_species]['total area'] += cell.volume
            plot_data[c_species]['count']      += 1      
        
        for cs in ['WT', 'MUT']:
            if plot_data[cs]['count'] != 0:
                plot_data[cs]['ave area'] = plot_data[c_species]['total area'] / plot_data[c_species]['count']
            
        return plot_data['WT']['ave area'], plot_data['MUT']['ave area']
        
        
        
    def step(self,mcs):
        global Extrusion_Data
        
        plot_data = {}
        for cs in ['WT', 'MUT']:
            plot_data[cs] = { 'ave area': 0. }
        plot_data['WT']['ave area'], plot_data['MUT']['ave area'] = self.get_wt_and_mut_average_area()            
            
        for cell in self.cell_list:                
            c_species = CELL_ATTR_TRACKER[str(cell.id)]['SPECIES']
            
            if cell.volume < (plot_data[c_species]['ave area'] * SIM_CONSTS[c_species]['extrusion threshold']):
                stress = (cell.volume - cell.targetVolume)**2.
                Extrusion_Data.append( [mcs, cell.id, c_species, cell.type, cell.volume, plot_data[c_species]['ave area'], stress] )
                self.delete_cell( cell )   
                    
        if SIM_CONSTS['SIM']['SAVE_DATA'] and (mcs % SIM_CONSTS['SIM']['SAVE_DATA_TO_CSV_INTERVAL']) and Extrusion_Data != []:
            update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Extrusion_Data_Filename, Extrusion_Data )   
            
            Extrusion_Data = []
                








 
class Plot_Total_Div_and_Apo_Counts(SteppableBasePy):        
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        global Tot_Div_Apo_Count_Data
        global Tot_Div_Apo_Count_Data_Filename
        Tot_Div_Apo_Count_Data_Filename = 'Total_Div_and_Apo_Counts.csv'        
        Tot_Div_Apo_Count_Col_Names = ['Time', 'WT Total Div Count', 'WT Total Apo Count', 'MUT Total Div Count', 'MUT Total Apo Count']
        
        if SIM_CONSTS['SIM']['SAVE_DATA']:
            file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Tot_Div_Apo_Count_Data_Filename, Tot_Div_Apo_Count_Col_Names)   
            if file_failed:
                kill_simulation(self, 'Failed to make ' + Tot_Div_Apo_Count_Data_Filename + 'Killing Sim') 
        Tot_Div_Apo_Count_Data = []
        
        
        if SIM_CONSTS['SIM']['SHOW_WT_DIV_APO_TOTALS_PLOT']:
            self.pW_wt_counts=self.add_new_plot_window( title = 'WT Division and Apoptosis Counts',
                                                 x_axis_title = 'Time (mcs)'           ,
                                                 y_axis_title = 'Count'          ,
                                                 x_scale_type = 'linear'               ,
                                                 y_scale_type = 'linear'               )       
            self.pW_wt_counts.add_plot('WT Div', style='Dots', color='green', size=5)   
            self.pW_wt_counts.add_plot('WT Apo', style='Dots', color='blue'  , size=5) 
            
        if SIM_CONSTS['SIM']['SHOW_MUT_DIV_APO_TOTALS_PLOT']:
            self.pW_mut_counts=self.add_new_plot_window( title = 'Mut Division and Apoptosis Counts',
                                                  x_axis_title = 'Time (mcs)'           ,
                                                  y_axis_title = 'Count'          ,
                                                  x_scale_type = 'linear'               ,
                                                  y_scale_type = 'linear'               )       
            self.pW_mut_counts.add_plot('MUT Div', style='Dots', color='purple' , size=5)   
            self.pW_mut_counts.add_plot('MUT Apo', style='Dots', color='orange' , size=3)            
            
            
            
    def step(self, mcs):
        global Tot_Div_Apo_Count_Data
        global TOTAL_DIV_AND_APO_COUNTS_editable
        

        plots_to_plot = []
        if SIM_CONSTS['SIM']['SHOW_WT_DIV_APO_TOTALS_PLOT']:
            plots_to_plot.append( [ 'WT', self.pW_wt_counts ] )            
        if SIM_CONSTS['SIM']['SHOW_MUT_DIV_APO_TOTALS_PLOT']:
            plots_to_plot.append( ['MUT', self.pW_mut_counts] )     

        #for cs,cs_plt in [['WT', self.pW_wt_counts], ['MUT', self.pW_mut_counts]]:
        if plots_to_plot != []:
            for cs,cs_plt in plots_to_plot:
                if SIM_CONSTS['SIM']['SHOW_'+cs+'_DIV_APO_TOTALS_PLOT']:
                    cs_plt.add_data_point( cs+' Div', mcs, TOTAL_DIV_AND_APO_COUNTS_editable[cs]['DIV'] )  
                    cs_plt.add_data_point( cs+' Apo', mcs, TOTAL_DIV_AND_APO_COUNTS_editable[cs]['APO'] )
                
            Tot_Div_Apo_Count_Data.append([ mcs, TOTAL_DIV_AND_APO_COUNTS_editable['WT']['DIV'],  TOTAL_DIV_AND_APO_COUNTS_editable['WT']['APO'], 
                                                TOTAL_DIV_AND_APO_COUNTS_editable['MUT']['DIV'], TOTAL_DIV_AND_APO_COUNTS_editable['MUT']['APO'] ])
                                                
            if SIM_CONSTS['SIM']['SAVE_DATA'] and (mcs % SIM_CONSTS['SIM']['SAVE_DATA_TO_CSV_INTERVAL']) and Tot_Div_Apo_Count_Data != []:
                update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Tot_Div_Apo_Count_Data_Filename, Tot_Div_Apo_Count_Data )
                Tot_Div_Apo_Count_Data = []
            


############
## Only use for pretty pictures (for science)
############
class Color_Cells_Switch_Between_GR_and_Cycle_Phase_Steppable(SteppableBasePy):
    def __init__(self,    frequency=1, 
                        color_GR_freq=20, color_Cell_Phase_mcs_offset=10,
                        kill_sim_after_colony_area_threshold_met = True,
                        kill_sim_mcs_buffer = 5000):
                        
        SteppableBasePy.__init__(self,   frequency)
        
        self.color_GR_freq               =  color_GR_freq
        self.color_Cell_Phase_mcs_offset =  color_Cell_Phase_mcs_offset
        
        self.kill_sim_after_colony_area_threshold_met = kill_sim_after_colony_area_threshold_met
        self.kill_sim_mcs_buffer                      = kill_sim_mcs_buffer
            
    def start(self):
        window_size = SIM_CONSTS['SIM']['X_MAX'] * SIM_CONSTS['SIM']['Y_MAX']
        self.colony_threshold_area = np.pi * ( ( SIM_CONSTS['SIM']['X_MAX']*0.5 ) - 25 )**2.
        
        self.colony_threshold_met  = False   
        self.colony_threshold_time = 100000000    

        self.check_for_colony_area_threshold = False
        self.check_for_colony_area_threshold = True
        
        
        
        global Colony_Threshold_size_Time_Data
        global Colony_Threshold_Filename
        Colony_Threshold_Filename = 'Colony_Size_Threshold_Time.csv'        
        Colony_Threshold_Col_Names = [ 'Time of Surpassing Colony Size Threshold']
        file_failed = make_new_csv_file(SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Colony_Threshold_Filename, Colony_Threshold_Col_Names)   
        
        if file_failed:
            kill_simulation(self, 'Failed to make ' + Colony_Threshold_Filename + 'Killing Sim')
        Colony_Threshold_size_Time_Data = []


    def set_cell_type_by_growth_rate(self):
        for cell in self.cell_list:

            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:

                dAdt = []
                for i,area_at_next_time in enumerate( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][1:] ):
                    dAdt.append( area_at_next_time - CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][i] )
                
                if len(dAdt) > 0:
                    #GAAAHHHHHHHH F#%$*% integer division!!  added 1.0 to fix
                    ave_G_rate = 1.0 * sum(dAdt) / len(dAdt)
                    
                    clr_scale = [ 0, 0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8 ]                    
                    no_change_in_type = True
                    
                    if ave_G_rate >= clr_scale[-1]:
                        cell.type = 11
                        no_change_in_type = False
                    elif ave_G_rate <= -clr_scale[-1]:
                        cell.type = 32
                        no_change_in_type = False
                    
                    if no_change_in_type:
                        if ave_G_rate <= 0:
                            abs_Gr = ave_G_rate * -1.
                        else:
                            abs_Gr = ave_G_rate * 1.
                        #abs_Gr = np.sqrt( ave_G_rate**2 )
                        
                        for j,s in enumerate( clr_scale[1:] ):
                        
                            if (clr_scale[j] <= abs_Gr) and (abs_Gr < s):
                                if ave_G_rate < 0:
                                    cell.type = 22 + j
                                    no_change_in_type = False
                                else:  # 
                                    cell.type = 21 - j
                                    no_change_in_type = False
    
    
    def update_recent_volumes_list(self):
        
        for cell in self.cell_list:
            if str(cell.id) not in CELL_ATTR_TRACKER['COLOR_CELLS']:
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] = []
        
        for cell in self.cell_list:                
            CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].append( cell.volume )
            
            if len( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] ) > 10:
                        
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].pop(0)  
                
        
        
    def set_cell_type_by_cycle_phase(self, mcs):
        for cell in self.cell_list:            
            
            if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
    
                if cell.volume <= CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]:
                    size_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]
                    
                    frac_G1_completed = float( abs(cell.volume - size_thresh*0.5) / (size_thresh*0.5) )
                    G1_phase_type = int( frac_G1_completed *4 ) 
                    
                    #cell.type = 20
                    if cell.type != (33 + G1_phase_type):
                        cell.type = 33 + G1_phase_type
                elif str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs']:
                    time_past = mcs - CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]
                    time_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str(cell.id)]
                    
                    frac_timer_phase_completed = 4 * float( time_past ) / float( time_thresh ) 
                    timer_phase_type = int( frac_timer_phase_completed )
                    
                    new_type = 37 + timer_phase_type                
                    if cell.type != new_type:
                        cell.type = new_type
                else:
                    cell.type = 37
            
        
        
        
    def step(self, mcs): 
        global Colony_Threshold_size_Time_Data
        
        self.update_recent_volumes_list()
        
        if mcs%self.color_GR_freq == 0:
            self.set_cell_type_by_growth_rate()
            
        if mcs%self.color_GR_freq == self.color_Cell_Phase_mcs_offset:
            self.set_cell_type_by_cycle_phase( mcs )
            
            
        
        
                
        
        
        if self.kill_sim_after_colony_area_threshold_met: 
            if (mcs%self.color_GR_freq == 0):
                total_area = 0
                for cell in self.cell_list:
                    total_area += cell.volume
                    
                if total_area >= self.colony_threshold_area:
                    self.check_for_colony_area_threshold = False
                    self.colony_threshold_time = mcs
                    self.colony_threshold_met  = True
                    
                    Colony_Threshold_size_Time_Data.append( [ mcs ] )
                    update_csv_file( SIM_CONSTS['SIM']['DATA_OUT_PATH'] + Colony_Threshold_Filename, Colony_Threshold_size_Time_Data )

            if self.colony_threshold_met:
                if mcs > (self.colony_threshold_time + self.kill_sim_mcs_buffer):
                    kill_simulation(self, 'Reached the Colony Size Threshold and past the Buffer time ...  Killing Sim')

                



############
## Only use for pretty pictures (for science)
############
class ColorCellsByGrowthRateSteppable(SteppableBasePy):
    def __init__(self,   frequency=1, update_color_freq=20):
        SteppableBasePy.__init__(self,   frequency)
        self.update_color_freq = update_color_freq
               
            
            
    def step(self, mcs):    
        
        for cell in self.cell_list:
            if str(cell.id) not in CELL_ATTR_TRACKER['COLOR_CELLS']:
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] = []
        
        for cell in self.cell_list:
            CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].append( cell.volume )
            
            #if CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE']:
            if len( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)] ) >= 10:
                dAdt = [] 
                if mcs%self.update_color_freq == 0:
                    if not( CELL_ATTR_TRACKER[str(cell.id)]['IS_ALIVE'] ):
                        cell.type = 9
                    else:
                        for i,area_at_next_time in enumerate( CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][1:] ):
                            dAdt.append( area_at_next_time - CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)][i] )
                        
                        #GAAAHHHHHHHH F#%$*% integer division!!  added 1.0 to fix
                        ave_G_rate = 1.0 *sum(dAdt) / len(dAdt)
                        
                        
                    
                    
                    
                        #clr_scale = [ 0, 0.25, 0.5, 1,    2, 3,   4, 5,   6, 7,   8 ]
                        #clr_scale = [ 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4,   5 ]
                        clr_scale = [ 0, 0.25, 0.5, 1,    2, 2.5, 3, 3.5, 4, 4.5, 5 ]
                        
                        no_change_in_type = True
                        
                        if ave_G_rate >= clr_scale[-1]:
                            cell.type = 11
                            no_change_in_type = False
                        elif ave_G_rate <= -clr_scale[-1]:
                            cell.type = 32
                            no_change_in_type = False
                        
                        if ave_G_rate <= 0:
                            abs_Gr = ave_G_rate * -1.
                        else:
                            abs_Gr = ave_G_rate * 1.
                        #abs_Gr = np.sqrt( ave_G_rate**2 )
                        
                        for i,s in enumerate( clr_scale[1:] ):
                        
                            if (clr_scale[i] <= abs_Gr) and (abs_Gr < s) and no_change_in_type:
                                if ave_G_rate < 0:
                                    cell.type = 22 + i
                                    no_change_in_type = False
                                else:  # 
                                    cell.type = 21 - i
                                    no_change_in_type = False
                    
                    
                        
                        
                #outside of the if mcs%10 --- delete oldest area value
                CELL_ATTR_TRACKER['COLOR_CELLS'][str(cell.id)].pop(0)    
                        
                        
                    
        
############
## Only use for pretty pictures (for science)
############ 
class ColorCellsByCellCyclePhaseSteppable(SteppableBasePy):
    def __init__(self,   frequency=1):
        SteppableBasePy.__init__(self,   frequency)
               
    def start(self):
        CELL_ATTR_TRACKER['COLOR_CELLS'] = {'Not_empty':[]}
            
    def step(self, mcs):    
        
        if mcs%5 == 0:
            for cell in self.cell_list:
                if cell.volume <= CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]:
                    size_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['sizer_thresholds'][str(cell.id)]
                    
                    frac_G1_completed = float( abs(cell.volume - size_thresh*0.5) / (size_thresh*0.5) )
                    G1_phase_type = int( frac_G1_completed *4 ) 
                    
                    #cell.type = 20
                    if cell.type != (33 + G1_phase_type):
                        cell.type = 33 + G1_phase_type
                elif str(cell.id) in CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs']:
                    time_past = mcs - CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_start_mcs'][str(cell.id)]
                    time_thresh = CELL_ATTR_TRACKER['MITOSIS_COND_DICT']['timer_thresholds'][str(cell.id)]
                    
                    frac_timer_phase_completed = 4 * float( time_past ) / float( time_thresh ) 
                    timer_phase_type = int( frac_timer_phase_completed )

                    if timer_phase_type > 4:
                        timer_phase_type = 4
                        kill_simulation(self, 'Error in ColorCellsByCellCyclePhaseSteppable() --- Killing Sim')
                    
                    new_type = 37 + timer_phase_type                
                    if cell.type != new_type:
                        cell.type = new_type
                else:
                    cell.type = 37
            
                        



class RadialForce(SteppableBasePy):
    def __init__(self,   frequency=1):
        SteppableBasePy.__init__(self,   frequency)
               
            
    def step(self, mcs):  

        for cell in self.cell_list:
            X, Y = abs(0.5*SIM_CONSTS['SIM']['X_MAX'] - cell.xCOM), abs(0.5*SIM_CONSTS['SIM']['Y_MAX'] - cell.yCOM)
            R = np.sqrt( X**2 + Y**2 )
            
            if (SIM_CONSTS['SIM']['EDGE_BUFFER'] < R) and (R < (0.5*SIM_CONSTS['SIM']['X_MAX'] - SIM_CONSTS['SIM']['EDGE_BUFFER'])):                               
                Fx_sign, Fy_sign = -1, -1
                if cell.xCOM < 0.5*SIM_CONSTS['SIM']['X_MAX']:
                    Fx_sign = 1
                if cell.yCOM < 0.5*SIM_CONSTS['SIM']['Y_MAX']:
                    Fy_sign = 1                    
                    
                ##F_r = F_MAGNITUDE * np.sin( (R - SIM_CONSTS['SIM']['EDGE_BUFFER']) * np.pi/(0.5*SIM_CONSTS['SIM']['X_MAX'] - 2*SIM_CONSTS['SIM']['EDGE_BUFFER']) )
                ##F_r = F_MAGNITUDE * np.cos( (R - SIM_CONSTS['SIM']['EDGE_BUFFER']) * (0.5*np.pi) /(0.5*SIM_CONSTS['SIM']['X_MAX'] - 2*SIM_CONSTS['SIM']['EDGE_BUFFER']) )
                F_r = F_MAGNITUDE
                                
                cell.lambdaVecX = Fx_sign * F_r * (X/R)
                cell.lambdaVecY = Fy_sign * F_r * (Y/R)
                
        # for cell in self.cell_list:
            # Fx_sign = 1.
            # if cell.xCOM > 0.5*SIM_CONSTS['SIM']['X_MAX']:
                # Fx_sign = -1.
            
            # cell.lambdaVecX = Fx_sign * F_MAGNITUDE




