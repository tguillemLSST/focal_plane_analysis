import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
import os

# Build GIF
#input_path = '/sps/lsst/users/tguillem/web/debug/images_dark_run/all3/13161/'
#input_path = '/sps/lsst/users/tguillem/web/batch/run5/dark/check_200522_flash/13161/'
input_path = '/sps/lsst/users/tguillem/web/batch/dark/ccd_assembly_fix_1/v4_1D/13161/after_master_bias/'
#000175/R14/S22/CCD_RawImage_2D_corr.png
os.makedirs(input_path+'gif', exist_ok = True) 
rafts=['R01' ,'R02' ,'R03' ,'R10' ,'R11' ,'R12' ,'R13' ,'R14' ,'R20' ,'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34' ,'R41' ,'R42' ,'R43']
#rafts_itl = ['R01' ,'R02' ,'R03' ,'R10' ,'R20', 'R41' ,'R42' ,'R43']
#rafts_e2v = ['R11' ,'R12' ,'R13' ,'R14', 'R21' ,'R22' ,'R23' ,'R24' ,'R30' ,'R31' ,'R32' ,'R33' ,'R34']
ccds=['S00' ,'S01' ,'S02' ,'S10' ,'S11' ,'S12' ,'S20' ,'S21' ,'S22']

for j in range(len(rafts)):
    for k in range(len(ccds)):
        with imageio.get_writer(input_path + 'gif/' + rafts[j] + '_' + ccds[k] + '.gif', mode='I', fps=1/2) as writer:
            for i in range(137,156):#155,160,165,170,175]:
                #filename = 'default'
                #if(i in range(5,10)):
                #    filename = input_path + '/bias_bias_00'+str(i)+'/R12/S22/CCD_Image_imutils_corr.png'
                #elif(i in range(20,39)):         
                #    filename = input_path + '/dark_dark_0'+str(i)+'/R12/S22/CCD_Image_imutils_corr.png'
                if(i==0):
                    filename = '/sps/lsst/users/tguillem/web/batch/overscan/dark/v4_1D/13161/after_master_bias/000140/' + rafts[j] + '/' + ccds[k] + '/CCD_RawImage_2D_corr.png'
                else:
                    filename = input_path + '/000'+str(i)+'/'+ rafts[j] + '/' + ccds[k] + '/CCD_RawImage_2D_corr.png'
                image = imageio.imread(filename)
                writer.append_data(image)        
        
