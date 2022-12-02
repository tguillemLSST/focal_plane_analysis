#!/usr/bin/env python
# coding: utf-8

##########
# Author: P. Antilogus
# Adapted by T. Guillemin 

# Goal: compute variance of super flats
##########

from eochar.bot_frame_op import *
from eochar.GetPhotoFlux import *

# system imports
import os
import time
import sys
from sys import exit
import glob

#  Specific package (display , pyfits ..)
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.visualization import ImageNormalize
from astropy.visualization import PercentileInterval
import matplotlib
matplotlib.rcParams['figure.dpi'] = 120

# color map "per amplifier"
cmap=plt.get_cmap('gist_ncar')
colors=[cmap(j)[:3] for j in np.linspace(0,1,17)]

# Data localisation
#BOT_REPO_DIR = '/sps/lsst/groups/FocalPlane/SLAC/LCA-10134_Cryostat-0001'
BOT_REPO_DIR = '/sps/lsst/groups/FocalPlane/SLAC/run5/'

#
# All existing rafts ======
#
raft_itl=['R01', 'R02', 'R03', 'R10', 'R20', 'R41', 'R42', 'R43']
raft_e2v=['R11', 'R12', 'R13', 'R14', 'R21', 'R22', 'R23', 'R24', 'R30', 'R31', 'R32', 'R33', 'R34']
raft_corner=['R00','R04','R40','R44']
# all existing sensor names
sensors_raft=['S00','S01','S02','S10','S11','S12','S20','S21','S22']
sensors_corner=['SG0','SG1','SW0','SW1']
sensors_8ch=['SW0','SW1']
# all channels  for sensors with 8 ( = sensors_8ch) or 16 channels
ch8=['C10','C11','C12','C13','C14','C15','C16','C17']
ch16=['C10','C11','C12','C13','C14','C15','C16','C17','C00','C01','C02','C03','C04','C05','C06','C07']
#
#
all_sensors={}
#
for raft in raft_itl+raft_e2v:
    all_sensors[raft]=sensors_raft
for raft in raft_corner:
    all_sensors[raft]=sensors_corner

# CONFIGURATION FOR THE CURRENT EXECUTION  ==========================================================================================================
# ---- raft and associated run ============ To be updated if needed 
# 'raft' :  list of raft 
# 'run' :  list of run , ex : ['9876','9874']  (remark : run is a string )   for all run of this raft : '*' 
#
run=['13144']
#raft=raft_corner+['R11','R14','R23','R33','R42']
raft=['R14']
#raft=['R14','R23','R33','R42']
#raft=['R14']
#raft=['R11','R14','R23','R33']
#raft=['R42']+raft_corner
#
# You can limit the list of sensor to consider to a sub list of the total raft sensor , or comment  the line bellow to keep them all 
all_sensors['R14']=['S22']
#
#directory to output data
output_data='/sps/lsst/users/tguillem/web/debug/'
#================= Method for the correction of bias 
# 1 Constant subtracted per image ( th constant is computed from the mean of the serial overscan ) 
#Bias_cor='Ct'
# 2D corection using individualy all line and serial overscan 
Bias_cor='2D' 
# directory patern with the image we want to use
# ================== Bias 
Flat=False
ImageLabel='Bias'
# Bias de réfrence
image_dir='flat_bias_5*'
# Bias entre flat 
#image_dir='flat_bias_*'
# flux max of the bis taken before the flat ...it's not perfect as the illumination is not uniform ...
flux_max=50000.
# ================== Flat
#Flat=True
#ImageLabel='Flat'
# flat serie
#image_dir='sflat_flat_SDSSi_H_*'
#
#covariance to be considered
# cov(0,0),cov(1,0),cov(0,1),cov(1,1) 
# ncovx=2
# ncovy=2
# Rahima , pour l'intant je propose que l'on oublie les convariance , on ne garde que la variance de la difference ============================
# cov(0,0)
ncovx=2
ncovy=2
# dont produce the table with the  overscan data (save space and time if not used ) 
over=False
# dont produce the error for the covariances (save space and time if not used ) 
cov_error=False
#cov_error=True
# ==================================================================================================================================================
ImageLabel+='_'+Bias_cor+'_bias_cor'
#
print("\n\nConfiguration summary ===============\n\n")
print("Run : ",run)
print("Raft: ",raft)
for raft_cur in raft : 
    print('Raft ',raft_cur,' with sensors ',all_sensors[raft_cur])
print("Plots will be written in : %s \n" % (output_data))
print("The bias correction will use the %s method. "% (Bias_cor))
if Flat : 
    print("The data read will be ploted  as  Flat with %s as label " % (ImageLabel))
else :
    print("The data read will be processed  as non-Flat with %s as label" % (ImageLabel))
    if image_dir=='flat_bias_*' :
        print ("The bias were taken after a pair of flat , we keep the BIAS only if the requested flux (from file name) of the previous flat is < %f " % (flux_max))
if over : 
    print('The overscan table will be filled ')
else : 
    print('The overscan table will not be filled ')
if cov_error:
    print('The error on covariance and image table will be filled ')
else : 
    print('The error on covariance and image table will not be filled ')
#    
print("We will compute covariances matrix up to cov(%d,%d) \n" % (ncovx,ncovy))    
print('The data are taken fron %s ' %(BOT_REPO_DIR ))
print('for the run(s) %s and from the directory %s\n\n' % (run,image_dir))

###function raft_param()
def raft_param(raft_cur):
        if raft_cur in raft_corner :
            nb_amp_in_raft=48
        else :
            nb_amp_in_raft=144
        #
        number_of_raft_sensor=len(all_sensors[raft_cur])
        # labels for the plots
        label_txt=np.zeros((number_of_raft_sensor),dtype=np.object_)
        label_pos=np.zeros((number_of_raft_sensor))
        if number_of_raft_sensor==9 : 
            label_chan=np.zeros((18),dtype=np.object_)
            label_chan_pos=np.zeros((18))
        else : 
            label_chan=np.zeros((6),dtype=np.object_)
            label_chan_pos=np.zeros((6))    
        ch_cur=0
        ichan_lab=0
        for  iccd in range(number_of_raft_sensor) : 
                label_chan[ichan_lab]='%d' % (1)
                label_chan_pos[ichan_lab]=ch_cur
                ichan_lab+=1
                if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                    number_of_amp_sensor=8
                else : 
                    number_of_amp_sensor=16
                    label_chan[ichan_lab]='%d' % (9)
                    label_chan_pos[ichan_lab]=ch_cur+8
                    ichan_lab+=1
                #
                label_txt[iccd]=all_sensors[raft_cur][iccd]
                label_pos[iccd]=ch_cur+number_of_amp_sensor/2-.5
                ch_cur+=number_of_amp_sensor
        return number_of_raft_sensor,nb_amp_in_raft,label_txt,label_pos,label_chan,label_chan_pos

###function SaveFig()
def SaveFig(fig,PlotFile,run_cur='',raft_cur='',ccd_cur='',hdu=0):
    if hdu<1 : 
        root_plt=os.path.join(output_data,run_cur,raft_cur,ccd_cur)
    else : 
        hdu_cur='%d' % (hdu)
        root_plt=os.path.join(output_data,run_cur,raft_cur,ccd_cur,hdu_cur)
    #        
    os.makedirs(root_plt,exist_ok=True)
    plotfile=os.path.join(root_plt,PlotFile)
    print ('PlotFile=',plotfile)
    fig.savefig(plotfile,bbox_inches='tight')
    plt.close(fig) 
    return

start_time=time.time()
###function PrintTimer()
def PrintTimer(text):
    sys.stdout.write("\r") 
    txt="%s ,  Running since %f seconds" % (text,time.time()-start_time)
    sys.stdout.write(txt) 
    sys.stdout.flush()

# loop on runs
for irun in range(len(run)) : 
    run_cur=run[irun]
    print('Analysis of run ',run_cur)
    # raw data directory for this run 
    # data=os.path.join(BOT_REPO_DIR,run_cur,'BOT_acq','v0')
    data=os.path.join(BOT_REPO_DIR,run_cur)
    print(data)
    all_raw=glob.glob(data+'/*')
    all_raw.sort()
    print(all_raw)
    # to select the last raw  data directory generated (in geneeral there is just 1 , bot some time more ...) 
    #data=all_raw[-1]
    print(data)
    # all image directories
    if image_dir=='flat_bias_*' :
        all_dir_0=glob.glob(data+'/'+image_dir)
        all_dir_0.sort()
        # look for flux in the flat taken just before 
        flat_dir=glob.glob(data+'/flat_*_flat1_*')
        flat1_data=[flat_dir[i].split('_') for i in range(len(flat_dir))]
        flat_bias_id=np.array([(int(flat1_data[i][-1])+1) for i in range(len(flat_dir))])
        flat_flux=np.array([(float(flat1_data[i][-3])) for i in range(len(flat_dir))])
        bias_arg_sort=np.argsort(flat_bias_id)
        flat_bias_id=flat_bias_id[bias_arg_sort]
        flat_flux=flat_flux[bias_arg_sort]
        #
        bias_id=np.array([int(all_dir_0[i].split('_')[-1]) for i in range(len(all_dir_0))])
        #
        pre_bias_flux=[]
        all_dir=[]
        # From there we can descide to keep only the bias taken after a given a flat with a flux bellow some value
        for i in range(0,len(all_dir_0)):
            if i==0 : 
                flux=0
            elif bias_id[i]!= flat_bias_id[i-1] :
                print ('ERROR bias_id %d # flat_bias_id %' % (bias_id[i],flat_bias_id[i-1]))
                break
            else :
                flux=flat_flux[i-1]
            # we selec bias taken after a low flux flat 
            if flux < flux_max :
                pre_bias_flux.append(flux)
                all_dir.append(all_dir_0[i])
        pre_bias_flux=np.array(pre_bias_flux)
    else :
        print('else')
        # all image directories
        all_dir=glob.glob(data+'/'+image_dir)
        all_dir.sort()
        print(all_dir)
        pre_bias_flux=np.zeros((len(all_dir)))
        #
    # fill the flux data
    nb_images=len(all_dir)
    #
    #
    for iraft in range(len(raft)) :
        raft_cur=raft[iraft]
        #
        #
        number_of_raft_sensor,nb_amp_in_raft,label_txt,label_pos,label_chan,label_chan_pos=raft_param(raft_cur)
        #
        #

        nb_pair_all=np.int32(nb_images/2)
        nb_files=nb_pair_all*2
        nb_amp=16
        flux=np.zeros((number_of_raft_sensor,nb_amp,nb_files))
        if over : 
            bias=np.zeros((number_of_raft_sensor,nb_files,nb_amp))
        dif_hist =np.zeros((number_of_raft_sensor,nb_amp,nb_files))
        dif_hist2=np.zeros((number_of_raft_sensor,nb_amp,nb_files))
        nb_pair=np.zeros((number_of_raft_sensor),dtype=np.int32)
        #
        for  iccd in range(number_of_raft_sensor) : 
            #
            if ( all_sensors[raft_cur][iccd] in sensors_8ch ):
                number_of_amp_sensor=8
            else : 
                number_of_amp_sensor=16
            #
            ccd_cur=all_sensors[raft_cur][iccd]
            #
            start=True
            #
            ifile=0
            i_pair=0
            file_cur=[]
            #
            for jfile in range(nb_images) :                
                if jfile==3:
                    break
                file_all=all_dir[jfile]+'/MC_C_*_'+raft_cur+'_'+ccd_cur+'.fits'
                file=glob.glob(file_all)
                print(file)
                if len(file)!=1 :
                    print( 'ERROR  : number of selected file incorect : ',file)
                    continue
                #ff=pyfits.open(file[0])
                file_cur.append(InFile(dirall=[file[0]],Slow=False,Bias=Bias_cor,verbose=False))
                if i_pair==0 : 
                    # check that the file is of the good LED 
                    #if not(int(file_cur[0][0].header['ILUM_LED']) in diode) :
                    #    file_cur[0].close()
                    #    file_cur=[]
                    #    nb_bad_photo+=1
                    #    continue
                    #print("1 single")
                    i_pair+=1
                    continue
                #
                if start : 
                    #
                    #print("1 =",tracemalloc.get_traced_memory()) 
                    xf=file_cur[0].all_file[0].first_col
                    xl=file_cur[0].all_file[0].first_s_over
                    yf=file_cur[0].all_file[0].first_line
                    yl=file_cur[0].all_file[0].first_p_over
                    nb_line=yl-yf
                    nb_col=xl-xf
                    print('+++++++++matrix')
                    print(xf)
                    print(xl)
                    print(yf)
                    print(yl)
                    print(file_cur[0].all_file[0])
                    #  
                    if over : 
                        overscan=np.zeros((nb_files,nb_amp),dtype=np.object_)
                    #
                    #
                    #
                    dif=np.zeros((nb_amp,nb_line,nb_col))
                    image=np.zeros((nb_amp,nb_line,nb_col))
                    prod=np.zeros((nb_amp,ncovy,ncovx,nb_line,nb_col))
                    if cov_error : 
                        image2=np.zeros((nb_amp,nb_line,nb_col))
                        prod2=np.zeros((nb_amp,ncovy,ncovx,nb_line,nb_col))
                        #
                        image2_offset=np.zeros((nb_amp))
                        prod2_offset=np.zeros((nb_amp,ncovy,ncovx))
                    #
                    start=False
                    print(str(nb_line)+' '+str(nb_col))
                    #2002 512
                #    
                # === 
                #print('*********DEBUG**********')
                for i in range(2) : 
                    for iamp in range(number_of_amp_sensor) : 
                        flux[iccd,iamp,ifile]=np.mean(np.median(file_cur[i].all_file[0].Image[iamp][yf+50:yl-20,xf+20:xl-20],axis=1),axis=0)
                        #print(file_cur[i].all_file[0].Image)
                        #print(np.shape(file_cur[i].all_file[0].Image))
                        #print(str(iamp)+' '+str(flux[iccd,iamp,ifile]))
                    ifile+=1
                    if ifile%10==0 : 
                        if Flat : 
                            PrintTimer(text='Run %s Raft %s CCD %s  %d Flats read '  % (run_cur,raft_cur,ccd_cur,ifile))
                        else : 
                            PrintTimer(text='Run %s Raft %s CCD %s  %d Bias read '  % (run_cur,raft_cur,ccd_cur,ifile))
                # compute per pair of images values
                for iamp in range(number_of_amp_sensor) : 
                    ima_1=(file_cur[0].all_file[0].Image[iamp][yf:yl,xf:xl])
                    ima_2=(file_cur[1].all_file[0].Image[iamp][yf:yl,xf:xl])
                    # 
                    # if Flats we have to egalize the flux ...if Bias we dont
                    if Flat : 
                        dif_cur=(ima_1/flux[iccd,iamp,ifile-2]-ima_2/flux[iccd,iamp,ifile-1])*(flux[iccd,iamp,ifile-2]+flux[iccd,iamp,ifile-1])/2.
                    else : 
                        dif_cur=(ima_1-ima_2)
                       
                    dif[iamp,:nb_line,:nb_col]+=dif_cur
                    # 
                    image[iamp,:nb_line,:nb_col]+=ima_1
                    image[iamp,:nb_line,:nb_col]+=ima_2
                    if cov_error : 
                        # offset the error, using the first pair ,  to avoid numerical problem 
                        if nb_pair[iccd]==0 : 
                            image2_offset[iamp]=np.mean(ima_1)**2
                        image2[iamp,:nb_line,:nb_col]+=ima_1**2-image2_offset[iamp]
                        image2[iamp,:nb_line,:nb_col]+=ima_2**2-image2_offset[iamp]
                    # 
                    dif_hist[iccd,iamp,nb_pair[iccd]]=dif_cur.mean()
                    dif_hist2[iccd,iamp,nb_pair[iccd]]=(dif_cur**2).mean()  
 
                    for iy in range(ncovy) :
                        for ix in range(ncovx) :
                            prod_cur=np.zeros((nb_line,nb_col))
                            if ix== 0 or iy== 0 : 
                                prod_cur[:nb_line-iy,:nb_col-ix]=dif_cur[:nb_line-iy,:nb_col-ix]*dif_cur[iy:nb_line,ix:nb_col]
                                prod[iamp,iy,ix,:,:]+=prod_cur
                                if cov_error : 
                                    if nb_pair[iccd]==0 : 
                                        prod2_offset[iamp,iy,ix]=np.mean(prod_cur)**2
                                    prod2[iamp,iy,ix,:,:]+=prod_cur**2-prod2_offset[iamp,iy,ix]
                                #del prod_cur
                            else :
                                prod2cur=np.zeros((nb_line,nb_col))
                                if ix*nb_line>iy*nb_col : 
                                    prod_cur[iy:nb_line-iy,:nb_col-ix]=dif_cur[iy:nb_line-iy,:nb_col-ix]*dif_cur[2*iy:nb_line,ix:nb_col]
                                    if cov_error : prod2cur[iy:nb_line-iy,:nb_col-ix]=dif_cur[iy:nb_line-iy,:nb_col-ix]*dif_cur[:nb_line-2*iy,ix:nb_col]
                                else :
                                    prod_cur[:nb_line-iy,ix:nb_col-ix]=dif_cur[:nb_line-iy,ix:nb_col-ix]*dif_cur[iy:nb_line,:nb_col-2*ix]
                                    if cov_error : prod2cur[:nb_line-iy,ix:nb_col-ix]=dif_cur[:nb_line-iy,ix:nb_col-ix]*dif_cur[iy:nb_line,2*ix:nb_col]
                                prod[iamp,iy,ix,:,:]+=(prod_cur+prod2cur)/2.
                                if cov_error : 
                                    # offset the error, using the first pair ,  to avoid numerical problem 
                                    if nb_pair[iccd]==0 : 
                                        prod2_offset[iamp,iy,ix]=np.mean(prod_cur)**2
                                    prod2[iamp,iy,ix,:,:]+=(prod_cur**2+prod2cur**2)/2.-prod2_offset[iamp,iy,ix]
                                #del prod_cur,prod2cur
                file_cur=[]
                i_pair=0
                nb_pair[iccd]+=1
            # average over file
            if nb_pair[iccd]==0 : continue 
            for iamp in range(number_of_amp_sensor) : 
                dif[iamp,:,:]/=nb_pair[iccd]
                image[iamp,:,:]/=nb_pair[iccd]*2
                if cov_error :
                    # compute the variance on image per pixels 
                    image2[iamp,:,:]=(image2[iamp,:,:]/(nb_pair[iccd]*2))+image2_offset[iamp]-image[iamp,:,:]**2   
                for iy in range(ncovy) :
                    for ix in range(ncovx) :
                        prod[iamp,iy,ix,:,:]/=nb_pair[iccd]
                        if cov_error : prod2[iamp,iy,ix,:,:]=(prod2[iamp,iy,ix,:,:]/nb_pair[iccd])-(prod[iamp,iy,ix,:,:]**2+prod2_offset[iamp,iy,ix])            
            
            # Plot per CCD ===============
            txt_common='for run %s raft %s ccd %s' %(run_cur,raft_cur,ccd_cur)
            # On dif images 
            fig=plt.figure(figsize=[15,20])
            txt='diff pair of %s: <axis=0>, val per column \n%s' % (ImageLabel,txt_common)
            plt.suptitle(txt)
            x=np.zeros((nb_amp))
            for i in range(nb_amp) :
                plt.subplot(16,1,i+1,title=i+1)
                x[i]=dif[i,50:-20,20:-20].mean()
                #plt.plot(dif[i,50:-20,20:-20].mean(axis=0))
                plt.plot(dif[i,:,:].mean(axis=0))
                plt.plot([0,nb_col],[0.,0.],color='g')
                txt='<>=%f' %(x[i])
                plt.plot([0,nb_col],[x[i],x[i]],color='r')
                plt.xlabel(txt)
            plt.show()
            SaveFig(fig,'dif_mean_axis0',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[15,20])
            txt='Diff pair of %s:<axis=1>, val per line \n%s' % (ImageLabel,txt_common)
            plt.suptitle(txt)
            for i in range(nb_amp) :
                plt.subplot(16,1,i+1,title=i+1)
                #plt.plot(dif[i,50:-20,20:-20].mean(axis=1))
                plt.plot(dif[i,:,:].mean(axis=1))
                plt.plot([0,nb_line],[0.,0.],color='g')
                txt='<>=%f' %(x[i])
                plt.plot([0,nb_line],[x[i],x[i]],color='r')
                plt.xlabel(txt)
            plt.show()
            SaveFig(fig,'dif_mean_axis1',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[15,15])
            title='Diff pair of %s: hist of pixels \n%s' % (ImageLabel,txt_common)
            for i in range(nb_amp) :
                plt.subplot(4,4,i+1,title=i+1)
                a=plt.hist(np.ravel(dif[i,100:-100,100:-100]),100)
            plt.show()
            #
            SaveFig(fig,'dif_hist',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[25,20])
            title='Stack diff pair of %s :  (70%s percentile) \n%s' % (ImageLabel,'%',txt_common)
            plt.suptitle(title)
            for i in range(nb_amp) :
                #norm = ImageNormalize(image[diode_cur,i,:,:]/x, interval=PercentileInterval(99.9),stretch=SqrtStretch())
                norm = ImageNormalize(dif[i,:,:], interval=PercentileInterval(70.))
                plt.subplot(2,8,i+1,title=i+1)
                #plt.imshow(dif[i,:,:],cmap = 'hot', vmin=-3, vmax=3, origin='lower',norm=norm)
                plt.imshow(dif[i,:,:],cmap = 'hot', vmin=-3, vmax=3, origin='lower')
                if not(i%8 ==0) :
                    figure=plt.gca()
                    y_axis = figure.axes.get_yaxis()
                    y_axis.set_visible(False)
                plt.colorbar()
            plt.show()
            SaveFig(fig,'dif_2d',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            #
            #fig=plt.figure(figsize=[25,20])
            #title='Variance for Stack diff pair of %s :  (70%s percentile) \n%s' % (ImageLabel,'%',txt_common)
            #plt.suptitle(title)
            #for i in range(nb_amp) :
            #    #norm = ImageNormalize(image[diode_cur,i,:,:]/x, interval=PercentileInterval(99.9),stretch=SqrtStretch())
            #    dif_var=prod[i,0,0,:,:]-(dif[i,:,:])**2
            #    norm = ImageNormalize(dif_var, interval=PercentileInterval(70.))
            #    plt.subplot(2,8,i+1,title=i+1)
            #    plt.imshow(dif_var,cmap = 'hot',origin='lower',norm=norm)
            #    if not(i%8 ==0) :
            #        figure=plt.gca()
            #        y_axis = figure.axes.get_yaxis()
            #        y_axis.set_visible(False)
            #    plt.colorbar()
            #plt.show()
            #del dif_var
            #SaveFig(fig,'dif_2d_var',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            
            # On images 
            fig=plt.figure(figsize=[15,20])
            txt='Super%s  uniformity: <axis=0>, val per column \n%s' % (ImageLabel,txt_common)
            plt.suptitle(txt)
            x=np.zeros((nb_amp))
            for i in range(nb_amp) :
                plt.subplot(16,1,i+1,title=i+1)
                x[i]=image[i,50:-20,20:-20].mean()
                txt='<>=%f' %(x[i])
                if Flat : 
                    #plt.plot(image[i,50:-20,20:-20].mean(axis=0)/x[i])
                    plt.plot(image[i,:,:].mean(axis=0)/x[i])
                    plt.plot([0,nb_col],[1.,1.],color='g')
                else : 
                    #plt.plot(image[i,50:-20,20:-20].mean(axis=0))
                    plt.plot(image[i,:,:].mean(axis=0))
                    plt.plot([0,nb_col],[x[i],x[i]],color='g')

                #plt.plot([0,nb_col],[x,x],color='r',label=txt)
                plt.xlabel(txt)
            plt.show()
            SaveFig(fig,ImageLabel+'_mean_axis0',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[15,20])
            txt='Super%s  uniformity:<axis=1>, val per line \n%s' % (ImageLabel,txt_common)
            plt.suptitle(txt)
            for i in range(nb_amp) :
                plt.subplot(16,1,i+1,title=i+1)
                txt='<>=%f' %(x[i])
                if Flat : 
                    #plt.plot(image[i,50:-20,20:-20].mean(axis=1)/x[i])
                    plt.plot(image[i,:,:].mean(axis=1)/x[i])
                    plt.plot([0,nb_line],[1.,1.],color='g')
                else : 
                    #plt.plot(image[i,50:-20,20:-20].mean(axis=1))
                    plt.plot(image[i,:,:].mean(axis=1))
                    plt.plot([0,nb_line],[x[i],x[i]],color='g')
                plt.xlabel(txt)
            plt.show()
            SaveFig(fig,ImageLabel+'_mean_axis1',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[15,15])
            title='Super%s  uniformity: hist of pixels \n%s' % (ImageLabel,txt_common)
            for i in range(nb_amp) :
                plt.subplot(4,4,i+1,title=i+1)
                a=plt.hist(np.ravel(image[i,100:-100,100:-100]),100)
            plt.show()
            #
            SaveFig(fig,ImageLabel+'_hist',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)
            fig=plt.figure(figsize=[25,20])
            title='Super%s uniformity (97%s percentile) \n %s' % (ImageLabel,'%',txt_common)
            plt.suptitle(title)
            for i in range(nb_amp) :
                #norm = ImageNormalize(image[diode_cur,i,:,:]/x, interval=PercentileInterval(99.9),stretch=SqrtStretch())
                norm = ImageNormalize(image[i,:,:], interval=PercentileInterval(97.))
                plt.subplot(2,8,i+1,title=i+1)
                plt.imshow(image[i,:,:],cmap = 'hot',origin='lower',norm=norm)
                if not(i%8 ==0) :
                    figure=plt.gca()
                    y_axis = figure.axes.get_yaxis()
                    y_axis.set_visible(False)
                plt.colorbar()
            plt.show()
            SaveFig(fig,ImageLabel+'_2d',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)  
            #
            fig=plt.figure(figsize=[25,20])
            title='Super%s variance uniformity (90%s percentile) \n %s' % (ImageLabel,'%',txt_common)
            plt.suptitle(title)
            for i in range(nb_amp) :
                #norm = ImageNormalize(image[diode_cur,i,:,:]/x, interval=PercentileInterval(99.9),stretch=SqrtStretch())
                norm = ImageNormalize(image2[i,:,:], interval=PercentileInterval(90.))
                plt.subplot(2,8,i+1,title=i+1)
                plt.imshow(image2[i,:,:],cmap = 'hot',origin='lower',norm=norm)
                if not(i%8 ==0) :
                    figure=plt.gca()
                    y_axis = figure.axes.get_yaxis()
                    y_axis.set_visible(False)
                plt.colorbar()
            plt.show()
            SaveFig(fig,ImageLabel+'_2d_var_90',run_cur=run_cur,raft_cur=raft_cur,ccd_cur=ccd_cur,hdu=0)            



