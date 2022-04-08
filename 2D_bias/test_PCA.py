import sys
import os
import shutil
import glob
import numpy as np
import lsst.eotest.sensor as sensorTest
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
from astropy.io import fits 

#job configuration
print('Configuration arguments: ', str(sys.argv))
str_run_all = str(sys.argv[1])
str_raft = str(sys.argv[2])
str_all_sensors = str(sys.argv[3])
repo_path=str(sys.argv[4])
root_dir=str(sys.argv[5])

#run_all=[str_run_all]
#raft=[str_raft]
#all_sensors[str_raft]=[str_all_sensors]
#butler = Butler(repo_path) 

#Run = '13159'
#det_name = 'R22_S22'
#bot_acq_dir = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/13162/BOT_acq/v0/107208'
#bot_acq_dir = '/sps/lsst/groups/FocalPlane/SLAC/run5/' + Run

Run = str_run_all
det_name = str_raft + '_' + str_all_sensors
bot_acq_dir = repo_path + Run
#hack for old run
#bot_acq_dir = repo_path

#flags
run_PCA = False

#bias_files = glob.glob(os.path.join(bot_acq_dir, 'flat_bias_5*', f'MC*{det_name}.fits'))
bias_files = glob.glob(os.path.join(bot_acq_dir, '*', f'MC*{det_name}.fits'))

#bias_files = [glob.glob(os.path.join(bot_acq_dir, 'bias_bias_012',f'MC*{det_name}.fits')),glob.glob(os.path.join(bot_acq_dir, 'bias_bias_013',f'MC*{det_name}.fits'))]
#bias_files = ['/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_012/MC_C_20211209_001682_R14_S12.fits', '/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_013/MC_C_20211209_001683_R14_S12.fits', '/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_014/MC_C_20211209_001684_R14_S12.fits']
#bias_files = ['/sps/lsst/groups/FocalPlane/SLAC/run5/13144/flat_bias_554/MC_C_20211207_000334_R14_S22.fits']
print(bias_files)
#output_path = '/sps/lsst/users/tguillem/web/PCA/superbias/'+Run+'/'+det_name+'/'
#output_path = '/sps/lsst/users/tguillem/web/PCA/flat/'+Run+'/'+det_name+'/'
output_path = root_dir + Run + '/' + det_name +'/'
if os.path.exists(output_path):
    shutil.rmtree(output_path)
os.makedirs(output_path)

# Prefix to use for output files.
file_prefix = root_dir + f'{det_name}_{Run}'

print('********************')
# Compute the PCA-based bias model, using bias_files as the input data.
if run_PCA==True:
    ccd_pcas = sensorTest.CCD_bias_PCA()
    pca_files = ccd_pcas.compute_pcas(bias_files, file_prefix)

#test overscan2D
#ccd_pcas = sensorTest.CCD_bias_overscan2D()
#file_prefix_local = f'{det_name}_{Run}'
#pca_files = ccd_pcas.compute_pcas(bias_files, file_prefix_local)

##########Check the PCA correction on bias in flat runs
#bias_files = ['/pbs/home/t/tguillem/data_run5/13144/flat_bias_200/MC_C_20211206_000755_R14_S12.fits']
#bias_files = ['/pbs/home/t/tguillem/data_run5/13144/flat_bias_749/MC_C_20211207_000529_R14_S12.fits', '/pbs/home/t/tguillem/data_run5/13144/flat_bias_515/MC_C_20211207_000295_R14_S12.fits']
#bias_files = glob.glob(os.path.join('/pbs/home/t/tguillem/data_run5/13144/', f'flat_bias_*',f'MC*{det_name}.fits'))
#print(bias_files)
#pca_files = ('R14_S22_13159_pca_bias.pickle', 'R14_S22_13159_pca_bias.fits')
#pca_files = ('R14_S12_13159_pca_bias.pickle', 'R14_S12_13159_pca_bias.fits')
#hack: mixed PCA files
pca_files = ('R14_S22_13159_pca_bias.pickle', 'R14_S22_13159_pca_bias.fits')
print('PCA-files used:' + str(pca_files))
##########

# Loop over bias files and amps in each frame and compute the
# bias-subtracted image using the PCA-bias model.  Print some stats
# for each residual image.
for i, bias_file in enumerate(bias_files):
    #if(i==5):
    #    break
    print(bias_file)
    frame = (bias_file.partition('/MC')[0])[-13:]
    #print(frame)

    print('+++++++PCA correction') 
    ccd = sensorTest.MaskedCCD(bias_file, bias_frame=pca_files)
    #print('ccd object OK')
    plt.figure(figsize=[25,20])
    plt.suptitle('PCA-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file)
    for amp in ccd:
        bias_subtracted_image = ccd.unbiased_and_trimmed_image(amp)
        imarr = bias_subtracted_image.getImage().array
        #print(i, amp, np.median(imarr),np.std(imarr))
        #plot
        #plt.figure()
        #plt.imshow(imarr.T, vmin=-10, vmax=10, cmap = 'hot', origin='lower')
        #plt.imshow(imarr, vmin=-3, vmax=3, cmap = 'hot', origin='lower')
        #plt.colorbar()
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        plt.colorbar()
        if not(amp%8==1) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
        #plt.imshow(image[i,:,:],cmap = 'hot',origin='lower',norm=norm)
        #if not(i%8 ==0) :
        #    figure=plt.gca()
        #    y_axis = figure.axes.get_yaxis()
        #    y_axis.set_visible(False)
        #   plt.colorbar()
    plt.savefig(output_path+"PCA_corr_bias_"+frame+".png") 
    plt.close()
    continue

    #images without PCA correction
    ccd_raw = sensorTest.MaskedCCD(bias_file)
    plt.figure(figsize=[25,20])
    plt.suptitle('Not corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file)
    for amp in ccd_raw:
        image = ccd_raw.unbiased_and_trimmed_image(amp)
        imarr = image.getImage().array
        #print(i, amp, np.median(imarr),np.std(imarr))
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr, vmin=-3, vmax=3, cmap = 'hot', origin='lower')
        plt.colorbar()
        if not(amp%8==1) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False) 
    plt.savefig(output_path+"bias_"+frame+".png")
    plt.close()
        
    print('+++++++2D correction')    
    #images with 2D raw correction
    #first_line,first_p_over,first_col,first_s_over=image_area(fitsfile)
    #self.first_col=first_col
    #self.first_s_over=first_s_over
    #self.first_line=first_line
    #self.first_p_over=first_p_over
    first_col=10
    first_s_over=522
    first_line=0
    first_p_over=2022
    last_l=2048
    last_s=576
    
    ccd_raw = sensorTest.MaskedCCD(bias_file)
    #HACK: get the master-bias corrected image from eotest turning off the PCA correction
    #ccd_raw = sensorTest.MaskedCCD(bias_file, bias_frame=pca_files)
    plt.figure(figsize=[25,20])
    plt.suptitle('2D-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file) 
    for amp in ccd_raw:
        image = ccd_raw.bias_subtracted_image(amp)
        imarr = image.getImage().array
        #print(imarr.shape)
        #########test 2D overscan correction
        # serial overscan
        mean_over_per_line=np.mean(imarr[:,first_s_over+2:],axis=1)
        #mean_over_per_line=np.zeros(2048)
        #####test prescan
        #for line in range(30,40):
            #print('line ' + str(line))
            #for col in range(0,10):
                #print(imarr[line][col])
        #mean_over_per_line_pre=np.mean(imarr[:,4:9],axis=1)
        #mean_over_per_line_post=np.mean(imarr[:,first_s_over+2:],axis=1)
        #mean_over_per_line_average=np.zeros(2048)
        #mean_over_sliding=np.zeros((2048,576))
        #for line in range(0,2047):
        #    mean_over_per_line_average[line]=
        #    #print(str(mean_over_per_line_pre[line])+' | ' + str(mean_over_per_line_post[line]))
        #    for col in range(0,575):
        #        mean_over_sliding[line][col]=(521-col)/512*mean_over_per_line_pre[line]+(col-10)/512*mean_over_per_line_post[line]
        #####
        #mean serial-parallel corner
        #mean_serial_parallel_corner = np.mean(imarr[first_p_over+2:,first_s_over+2:])
        #print(mean_serial_parallel_corner)
        #print('format mean_over_per_line ' + str(mean_over_per_line.shape))
        rawl=np.zeros((last_l-first_p_over-2,last_s))
        # // ovesrcan per column , corrected by the serial value per line
        for l in range(first_p_over+2,last_l):
            rawl[l-first_p_over-2,:]=imarr[l,:]-mean_over_per_line[l]
            #rawl[l-first_p_over-2,:]=imarr[l,:]-mean_serial_parallel_corner
        #print('format rawl ' + str(rawl.shape))     
        # // overscan
        mean_over_per_column=np.mean(rawl[:,:],axis=0)
        #print('format mean_over_per_column ' + str(mean_over_per_column.shape))
        # // correction
        linef=mean_over_per_line[:,np.newaxis]
        #print('format linef ' + str(linef.shape))
        # generate the 2D correction (thank's to numpy)
        over_cor_mean=mean_over_per_column+linef
        #print('format over_cor_mean ' + str(over_cor_mean.shape))
        # 2D correction of the overscan : 1 overscan subtracted per line , 1 overscan subtracted per column
        imarr -= over_cor_mean
        #print(over_cor_mean[100][100])
        #print(linef[100])
        #print(mean_over_per_column[100])
        #print('===========zero')
        #imarr -= mean_over_sliding
        
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        plt.colorbar()
        if not(amp%8==1) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
            #plt.imshow(image[i,:,:],cmap = 'hot',origin='lower',norm=norm)
            #if not(i%8 ==0) :
            #    figure=plt.gca()
            #    y_axis = figure.axes.get_yaxis()
            #    y_axis.set_visible(False)
            #   plt.colorbar()
    plt.savefig(output_path+"2D_corr_bias_"+frame+".png") 
    plt.close()

    #2D sliding window
    #print('+++++++2D sliding window correction') 
    #plt.suptitle('2D-sliding-window-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file) 
    #for amp in ccd_raw:
    #    image = ccd_raw.bias_subtracted_image(amp)
    #    imarr = image.getImage().array
    #    #print(imarr.shape)
    #    #########test 2D overscan correction
    #    # serial overscan
    #    mean_over_per_line=np.mean(imarr[:,first_s_over+2:],axis=1)
    #    #print('format mean_over_per_line ' + str(mean_over_per_line.shape))
    #    rawl=np.zeros((last_l-first_p_over-2,last_s))
    #    # // ovesrcan per column , corrected by the serial value per line
    #    for l in range(first_p_over+2,last_l):
    #        rawl[l-first_p_over-2,:]=imarr[l,:]-mean_over_per_line[l]
    #    #print('format rawl ' + str(rawl.shape))     
    #    # // overscan
    #    mean_over_per_column=np.mean(rawl[:,:],axis=0)
    #    #print('format mean_over_per_column ' + str(mean_over_per_column.shape))
    #    # // correction
    #    linef=mean_over_per_line[:,np.newaxis]
    #    #print('format linef ' + str(linef.shape))
    #    # generate the 2D correction (thank's to numpy)
    #    over_cor_mean=mean_over_per_column+linef
    #    #print('format over_cor_mean ' + str(over_cor_mean.shape))
    #    # 2D correction of the overscan : 1 overscan subtracted per line , 1 overscan subtracted per column
    #    mean_over_sliding=np.zeros((2048,576))
    #    for line in range(0,2047):
    #        for col in range(0,575):
    #            mean_over_sliding[line][col]=np.mean(imarr[0:100,first_s_over+2:])
    #    
    #    imarr -= mean_over_sliding
    #    #print(over_cor_mean[100][100])
    #    #print(linef[100])
    #    #print(mean_over_per_column[100])
    #    #print('===========')
    #    
    #    plt.subplot(2,8,amp,title=amp)
    #    plt.imshow(imarr, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
    #    plt.colorbar()
    #    if not(amp%8==1) :
    #        figure=plt.gca()
    #       y_axis = figure.axes.get_yaxis()
    #        y_axis.set_visible(False)
    #        #plt.imshow(image[i,:,:],cmap = 'hot',origin='lower',norm=norm)
    #        #if not(i%8 ==0) :
    #        #    figure=plt.gca()
    #        #    y_axis = figure.axes.get_yaxis()
    #        #    y_axis.set_visible(False)
    #        #   plt.colorbar()
    #plt.savefig(output_path+"2D_sliding_corr_bias_"+frame+".png") 
    #plt.close()
    
    #line
    print('+++++++1D line correction') 
    plt.figure(figsize=[25,20])
    plt.suptitle('1D-line-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file)
    for amp in ccd_raw:
        image = ccd_raw.bias_subtracted_image(amp)
        imarr_line = image.getImage().array

        #1D line overscan
        mean_over_per_line=np.mean(imarr_line[:,first_s_over+2:],axis=1)
 
        line_2D = np.zeros((2048,576))
        for col in range(0,575):
            line_2D[:,col]=mean_over_per_line[:]
        #print(line_2D)
        imarr_line -= line_2D
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr_line, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        plt.colorbar()
        if not(amp%8==1) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
            #plt.imshow(image[i,:,:],cmap = 'hot',origin='lower',norm=norm)
            #if not(i%8 ==0) :
            #    figure=plt.gca()
            #    y_axis = figure.axes.get_yaxis()
            #    y_axis.set_visible(False)
            #   plt.colorbar()
    plt.savefig(output_path+"1D_line_corr_bias_"+frame+".png") 
    plt.close()
    
    #column
    print('+++++++1D column correction') 
    plt.figure(figsize=[25,20])
    plt.suptitle('1D-column-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file)
    for amp in ccd_raw:
        image = ccd_raw.bias_subtracted_image(amp)
        imarr_col = image.getImage().array

        #repeated, to improve
        rawl=np.zeros((last_l-first_p_over-2,last_s))
        for l in range(first_p_over+2,last_l):
            rawl[l-first_p_over-2,:]=imarr_col[l,:]
        mean_over_per_column=np.mean(rawl[:,:],axis=0)     
        
        col_2D = np.zeros((2048,576))
        for line in range(0,2047):
            col_2D[line,:]=mean_over_per_column[:]
        imarr_col -= col_2D
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr_col, vmin=-5, vmax=5, cmap = 'hot', origin='lower')
        plt.colorbar()
        if not(amp%8==1) :
            figure=plt.gca()
            y_axis = figure.axes.get_yaxis()
            y_axis.set_visible(False)
    plt.savefig(output_path+"1D_column_bias_"+frame+".png") 
    plt.close()
    
    #plot mean_overscan_column vs column number

    
#check pca_superbias function
#superbias = sensorTest.pca_superbias(bias_files, pca_files, outfile='superbias.FITS', overwrite=True,statistic=afwMath.MEDIAN)

print('DONE')
