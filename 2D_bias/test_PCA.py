import os
import glob
import numpy as np
import lsst.eotest.sensor as sensorTest
import matplotlib.pyplot as plt

Run = '13145'
det_name = 'R14_S12'
#bot_acq_dir = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/13162/BOT_acq/v0/107208'
bot_acq_dir = '/sps/lsst/groups/FocalPlane/SLAC/run5/13151/'

bias_files = glob.glob(os.path.join(bot_acq_dir, 'bias_bias*',
                                    f'MC*{det_name}.fits'))
#bias_files = [glob.glob(os.path.join(bot_acq_dir, 'bias_bias_012',f'MC*{det_name}.fits')),glob.glob(os.path.join(bot_acq_dir, 'bias_bias_013',f'MC*{det_name}.fits'))]
#bias_files = ['/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_012/MC_C_20211209_001682_R14_S12.fits', '/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_013/MC_C_20211209_001683_R14_S12.fits', '/sps/lsst/groups/FocalPlane/SLAC/run5/13151/bias_bias_014/MC_C_20211209_001684_R14_S12.fits']
print(bias_files)
output_path = '/sps/lsst/users/tguillem/web/PCA/'+Run+'/'+det_name+'/'

# Prefix to use for output files.
file_prefix = f'{det_name}_{Run}'

print('********************')
# Compute the PCA-based bias model, using bias_files as the input data.
ccd_pcas = sensorTest.CCD_bias_PCA()
pca_files = ccd_pcas.compute_pcas(bias_files, file_prefix)
print(pca_files)

# Loop over bias files and amps in each frame and compute the
# bias-subtracted image using the PCA-bias model.  Print some stats
# for each residual image.
for i, bias_file in enumerate(bias_files):
    if(i==2):
        break
    print(bias_file)
    ccd = sensorTest.MaskedCCD(bias_file, bias_frame=pca_files)
    plt.figure(figsize=[25,20])
    plt.suptitle('PCA-corrected image: run '+ Run + ' detector ' + det_name + '\n for file ' + bias_file)
    for amp in ccd:
        bias_subtracted_image = ccd.unbiased_and_trimmed_image(amp)
        imarr = bias_subtracted_image.getImage().array
        print(i, amp, np.median(imarr),np.std(imarr))
        #plot
        #plt.figure()
        #plt.imshow(imarr.T, vmin=-10, vmax=10, cmap = 'hot', origin='lower')
        #plt.imshow(imarr, vmin=-3, vmax=3, cmap = 'hot', origin='lower')
        #plt.colorbar()
        plt.subplot(2,8,amp,title=amp)
        plt.imshow(imarr, vmin=-3, vmax=3, cmap = 'hot', origin='lower')
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
    plt.savefig(output_path+"PCA_corr_bias_"+str(i)+".png") 
    
