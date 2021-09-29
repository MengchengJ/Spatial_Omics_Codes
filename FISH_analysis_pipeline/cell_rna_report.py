import os
import numpy as np
import pandas as pd
import cv2
import yaml
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "pgf.texsystem": "xelatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'figure.dpi': 300,
})
from skimage.feature import peak_local_max
from skimage import img_as_float
from scipy.io import loadmat

ellipse_kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(7,7))

def tophat_spots(image, kernel):
    return cv2.morphologyEx(image, cv2.MORPH_TOPHAT, kernel)

def extract_coordinates(image, kernel, savename, snr=4, quantile=0.25, smooth=False):
    if smooth:
        image = cv2.blur(image,(3,3))
    image_tophat = img_as_float(tophat_spots(image,kernel))
    coordinates = peak_local_max(image_tophat)
    coordinates = coordinates[image[coordinates[:,0],coordinates[:,1]]>snr*np.mean(image)]
    intensities = image[coordinates[:,0],coordinates[:,1]]
    plt.hist(intensities,bins=1000)
    plt.xlabel('Intensity')
    plt.ylabel('Number of cells')
    plt.title('Intensity Distribution of Spots')
    plt.savefig(savename)
    plt.close()
    threshold = np.quantile(intensities,quantile)
    coordinates = coordinates[image[coordinates[:,0],coordinates[:,1]]>threshold]
    return coordinates,threshold

def nearest_center(coordinates,cell_pos,im):
    assigned_label = []
    for i in range(coordinates.shape[0]):
        if np.min(np.linalg.norm(cell_pos-coordinates[i],axis=1)) <= 50:
            assigned_label.append(np.argmin(np.linalg.norm(cell_pos-coordinates[i],axis=1)))
        else:
            assigned_label.append(-1)
    intensity_table = pd.DataFrame({'X':coordinates[:,1],'Y':coordinates[:,0],'intensity':im[coordinates[:,0],coordinates[:,1]],'label':assigned_label})
    return intensity_table

def addlabels(x,y):
    bleed = 0.01 * y.max()
    for i in range(len(x)):
        plt.text(x[i], y.iloc[i]+bleed, y.iloc[i], ha = 'center')

def bar_rna_count(int_table,_label,add_label=False,start=0):
    label_count = int_table.groupby('label').count()
    unassigned = int(label_count.iloc[0,0])
    label_count = label_count.loc[1:,:]
    cell_counts = label_count.groupby('intensity').count()
    cell_counts = cell_counts[cell_counts.index <= 20]
    plt.bar(cell_counts.index[start:],cell_counts['X'][start:],alpha=1,label=_label)
    if add_label:
        addlabels(cell_counts.index[start:],cell_counts['X'][start:])
    return unassigned

BASE_DEST_DIRECTORY = r'D:\FISH_images_processed'
def generate_report(run_id,cycle,channel,snr):
    report_dict = {}
    run_id_p = f'{run_id}_processed'
    label_img_mat = loadmat(os.path.join(BASE_DEST_DIRECTORY,run_id_p,'CellMap.mat'))
    label_img = label_img_mat['CellMap']
    cell_pos_mat = loadmat(os.path.join(BASE_DEST_DIRECTORY,run_id_p,'CellYX.mat'))
    cell_pos = cell_pos_mat['CellYX']
    im = cv2.imread(os.path.join(BASE_DEST_DIRECTORY,run_id_p,'stitched',f'cyc_{cycle}_{channel}.tif'), -cv2.IMREAD_ANYDEPTH)
    report_dict['Image mean intensity'] = float(np.mean(im))
    coordinates,threshold = extract_coordinates(im,ellipse_kernel,os.path.join(BASE_DEST_DIRECTORY,run_id_p,'report',f'{run_id}_cyc_{cycle}_{channel}_int_dist.pdf'),snr)
    report_dict['Spot count'] = int(coordinates.shape[0])
    report_dict['Threshold'] = float(threshold)
    int_table_cell_map = pd.DataFrame({'X':coordinates[:,1],'Y':coordinates[:,0],'intensity':im[coordinates[:,0],coordinates[:,1]],'label':label_img[coordinates[:,0],coordinates[:,1]]})
    int_table_nearest_center = nearest_center(coordinates,cell_pos,im)
    unassigned_nearest_center = bar_rna_count(int_table_nearest_center,'Nearest center',add_label=True)
    unassigned_cell_map = bar_rna_count(int_table_cell_map,'Cell map')
    report_dict['Unassigned spot count (nearest center)'] = unassigned_nearest_center
    report_dict['Unassigned spot count (cell map)'] = unassigned_cell_map
    plt.xlabel('Cell RNA counts')
    plt.ylabel('Number of cells')
    plt.title('Distribution of Spots in Individual Cell')
    plt.legend()
    plt.savefig(os.path.join(BASE_DEST_DIRECTORY,run_id_p,'report',f'{run_id}_cyc_{cycle}_{channel}_spot_cell_count.pdf'))
    plt.close()
    with open(os.path.join(BASE_DEST_DIRECTORY,run_id_p,'report',f'{run_id}_cyc_{cycle}_{channel}_report.yml'),'w') as f:
        yaml.dump(report_dict, f, default_flow_style=False)

if __name__ == "__main__":
    generate_report('20210722-1-ShortenRCAtest',2,'cy3',2)