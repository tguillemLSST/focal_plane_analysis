import os
import glob
import matplotlib.pyplot as plt
from lsst.eotest.sensor import CCD_bias_PCA
plt.ion()

#ccd = 'R22_S22'
#amp = 3
ccd = 'R14_S21'
amp = 8
#pickle_files = sorted(glob.glob(f'{ccd}_amp{amp:02d}_*_30_30_*10_4_pca_bias.pickle'))
#file = open('/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/eotest/python/lsst/eotest/sensor/R22_S11_13151_pca_bias.pickle',"rb" )
#pickle_files = sorted(glob.glob('/sps/lsst/users/tguillem/Rubin/Focal_Plane/lsst_distrib/w_2022_01/eotest/python/lsst/eotest/sensor/R22_S11_13151_pca_bias.pickle'))
pickle_files = ['/pbs/home/t/tguillem/web/PCA/batch/R14_S22_13159_pca_bias.pickle']
#pickle_files = ['R14_S22_13159_pca_bias.pickle']

pickle_title = 'R14_S22_13159_pca_bias'
for pickle_file in pickle_files:
    print(pickle_file)
    #ncomps = [int(_) for _ in pickle_file.split('_')[4:6]]
    ccd_pcas = CCD_bias_PCA.read_pickle(pickle_file)
    pcax, pcay = ccd_pcas[amp]
    pca_dict = dict(pcax=pcax, pcay=pcay)
    for pca_label, pca in pca_dict.items():
        ncomp = pca.n_components
        plt.figure(figsize=(9, 4*(ncomp//3)))
        for i in range(ncomp):
            #variance
            plt.subplot(4, 2, i+1)
            plt.plot(pca.components_[i])
            explained_variance = pca.explained_variance_ratio_[i]*100
            plt.title(f'var: {explained_variance:.3f}%')
            #            plt.title(f'amp {amp}, {pca_label} component {i+1}, '
            #                      f'explained variance: {explained_variance:.3f}%')
            #plt.tight_layout(rect=(0, 0, 1, 0.97))
            #plt.ylim([-0.1, 0.1])
            plt.suptitle(pickle_title)
            #file_prefix = pickle_file[:len(f'{ccd}_amp{amp:02d}_12845_06_08')]
            #plt.savefig(f'{pca_label}_{file_prefix}_PCA_components.png')
            plt.subplots_adjust(left=0.1,
                                bottom=0.1,
                                right=0.9,
                                top=0.9,
                                wspace=0.4,
                                hspace=0.4)
            plt.savefig(f'{pca_label}_components.png')

            #variance
            #plt.subplot(ncomp, 3, i+1)
            #plt.plot(pca.explained_variance_[i])
            #explained_variance = pca.explained_variance_ratio_[i]*100
            ##            plt.title(f'amp {amp}, {pca_label} component {i+1}, '
            ##                      f'explained variance: {explained_variance:.3f}%')
            #plt.tight_layout(rect=(0, 0, 1, 0.97))
            #plt.suptitle(pickle_title)
            ##file_prefix = pickle_file[:len(f'{ccd}_amp{amp:02d}_12845_06_08')]
            ##plt.savefig(f'{pca_label}_{file_prefix}_PCA_components.png')
            #plt.savefig(f'Variance_{pca_label}_components.png')
