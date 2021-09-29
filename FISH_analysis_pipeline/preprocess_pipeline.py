import os
from os_snippets import remove_empty_folders
import numpy as np
import pandas as pd
import cv2
from scipy.io import loadmat
from skimage.io import imread
from skimage.io import imsave
from skimage.util import img_as_uint
from glob import glob
from tqdm import tqdm
import matlab.engine
from all_in_focus import all_in_focus
from focal_stack_legacy import focal_stack
from dft_registration import register_with_ref
from dft_registration import register_channel
from dft_registration import interchannel_correction
from dft_registration import register_pair
from MIST_run_ijm import create_mist_command
from MIST_run_ijm import run_ijm
from stitch_meta import stitch_from_meta
from cell_segmentation import segment_cell
from cell_rna_report import generate_report
from report_xetex import integrate_report

BASE_DIRECTORY = r'D:\FISH_images'
BASE_DEST_DIRECTORY = r'E:\FISH_images_processed'
FSTACK_DIRECTORY = 'C:/Users/Dell/Documents/LabView_FISH_PKU/FISH_analysis_pipeline/fstack'
CIDRE_DIRECTORY = 'C:/Users/Dell/Documents/LabView_FISH_PKU/FISH_analysis_pipeline/CIDRE'
CHANNELS = ['cy3','cy5','FAM','DAPI']

def try_mkdir(d):
    """Try to make a new directory d. Passes if it already exists."""
    try:
        os.makedirs(d)
    except FileExistsError:
        print(f'{d} exists.')
        pass

def imlist_to_df(imlist):
    """Given a list of properly named .tif files, return a dataframe with organized cycle, tile, channel and depth information."""
    imlist_2d = [s.strip('.tif').split('-') for s in imlist]
    df = pd.DataFrame(imlist_2d, columns=['Cycle','Tile','Channel','Z'])
    df = df.loc[:,['Cycle','Tile','Z']]
    for column in ['Cycle','Tile','Z']:
        df[column] = df[column].apply(lambda x: int(x[1:]))
    return df

def get_stack_num(path):
    """Given a directory, return max z-axis index plus one (which is the image count of a stack)."""
    sample = [f for f in os.listdir(path) if f.endswith('.tif')][0]
    prefix = sample.split('-Z')[0]
    stack_num = max([int(f.strip('.tif').split('-Z')[1]) for f in os.listdir(path) if f.startswith(prefix)]) + 1
    return stack_num

def get_tiles(d):
    """Given a directory with tile information, return a series of tile index."""
    try:
        tile_info = loadmat(os.path.join(d,'TileInfo.mat'))
    except FileNotFoundError as exc:
        raise FileNotFoundError('Tile info file not included') from exc
    tile_x = int(tile_info['TileX'])
    tile_y = int(tile_info['TileY'])
    tile_num = tile_x * tile_y
    tiles = [i+1 for i in range(tile_num)]
    return tile_x,tile_y,tiles

def get_cycles(d,even_less=False,zero=False):
    """Given a directory, return full path of folders named cyc_*.
    
    Keyword arguments:
    even_less -- only odd cycles will be included (default False)
    zero -- include cycle 0 (default False)
    """
    if not even_less:
        cycles = glob(os.path.join(d,'cyc_*'))
        cycles = [c for c in cycles if not c.endswith('cyc_0')]
    else:
        cycles = glob(os.path.join(d,'cyc_*[13579]'))
        if zero:
            cycles += glob(os.path.join(d,'cyc_0'))
    #print(cycles)
    return cycles

def image_stacking(in_directory,out_directory,even_less=False,zero=False):
    """Full all-in-focus image generation pipeline, given Run ID."""
    try_mkdir(out_directory)
    cycles = get_cycles(in_directory,even_less,zero)
    stack_num = get_stack_num(cycles[0])
    chn_info = {}
    for chn in CHANNELS:
        chn_info[chn] = imlist_to_df([f for f in os.listdir(cycles[0]) if chn in f])
        if glob(os.path.join(in_directory,'cyc_*',f'*-{chn}-*')):
            eng = matlab.engine.start_matlab()
            eng.cd(FSTACK_DIRECTORY)
            eng.fstack_stacknum_channel(in_directory,out_directory,stack_num,chn,nargout=0)
            eng.quit()
    remove_empty_folders(out_directory)

def shading_correction(in_directory,out_directory):
    """Correct background shading of each child directory individually, given parent directory."""
    try_mkdir(out_directory)
    eng = matlab.engine.start_matlab()
    eng.cd(CIDRE_DIRECTORY)
    eng.background_correction(in_directory,out_directory,nargout=0)
    eng.quit()

def registration_batch(in_directory,out_directory,src_directory,ref_chn='cy3',ref_cyc=1):
    """Register following cycles using a particular channel & cycle as reference."""
    try_mkdir(out_directory)
    chn_info = {}
    for chn in CHANNELS:
        chn_info[chn] = imlist_to_df([f for f in os.listdir(get_cycles(src_directory)[0]) if chn in f and not f.startswith('.')])
    captured_tiles = list(chn_info[ref_chn]['Tile'].unique())
    img_names = [f'FocalStack_{t:03d}.tif' for t in captured_tiles]
    register_with_ref(in_directory,out_directory,CHANNELS,img_names,ref_cyc,ref_chn)

def registration_separate(in_directory,out_directory,src_directory,ref_chn='cy3',ref_cyc=1):
    """Correct shifts between multiple channels, then align each channel separately."""
    try_mkdir(out_directory)
    chn_info = {}
    for chn in CHANNELS:
        chn_info[chn] = imlist_to_df([f for f in os.listdir(get_cycles(src_directory)[0]) if chn in f and not f.startswith('.')])
    captured_tiles = list(chn_info[ref_chn]['Tile'].unique())
    img_names = [f'FocalStack_{t:03d}.tif' for t in captured_tiles]
    interchannel_correction(in_directory,out_directory,ref_cyc,ref_chn,img_names)
    for chn in CHANNELS:
        if chn in chn_info:
            register_channel(in_directory,out_directory,ref_cyc,chn,img_names)

def patch_missing_tile(d,tiles):
    """Patch (possible) missing images in a given tile grid."""
    cyc_chn_list = [f for f in os.listdir(d) if f.startswith('cyc_')]
    # print(glob(os.path.join(d,cyc_chn_list[0],'*.tif')))
    sample = imread(glob(os.path.join(d,cyc_chn_list[0],'*.tif'))[0])
    empty_im = np.zeros(sample.shape,dtype=np.uint16)
    for cyc_chn in cyc_chn_list:
        imgs = [f for f in os.listdir(os.path.join(d,cyc_chn)) if f.endswith('.tif')]
        for tile in tiles:
            name = f'FocalStack_{tile:03d}.tif'
            if name not in imgs:
                imsave(os.path.join(d,cyc_chn,name),empty_im,check_contrast=False)

def template_stitch(in_directory,out_directory,tile_x,tile_y,ref_chn='cy3',ref_cyc=1):
    """MIST stitching of images in template channel & cycle (unregistered)."""
    try_mkdir(out_directory)
    template_cyc_chn = f'cyc_{ref_cyc}_{ref_chn}'
    create_mist_command(tile_x,tile_y,os.path.join(in_directory,template_cyc_chn).replace('\\',r'\\\\'),out_directory.replace('\\',r'\\\\'))
    if glob(os.path.join(out_directory,f'*global*')):
        pass
    else:
        run_ijm('MIST_temp.ijm')

def segment_save(in_directory,out_directory,cyc='10',chn='DAPI'):
    """Try segment stitched image given cycle & channel, then save the results."""
    cell_im = imread(os.path.join(in_directory,f'cyc_{cyc}_{chn}.tif'))
    coordinates,segmented = segment_cell(cell_im)
    np.savetxt(os.path.join(out_directory,'cell_coordinates.txt'),coordinates,fmt='%d')
    imsave(os.path.join(out_directory,'cell_segmented.tif'),segmented,check_contrast=False)

def main(run_id,segment=False):
    """Full pipeline of image preprocessing."""
    raw_directory = os.path.join(BASE_DIRECTORY,run_id)
    dest_directory = os.path.join(BASE_DEST_DIRECTORY,f'{run_id}_processed')
    try_mkdir(dest_directory)
    aif_directory = os.path.join(dest_directory,'focal_stacked')
    sdc_directory = os.path.join(dest_directory,'background_corrected')
    rgs_directory = os.path.join(dest_directory,'registered')
    stc_directory = os.path.join(dest_directory,'stitched')
    report_directory = os.path.join(dest_directory,'report')
    image_stacking(raw_directory,aif_directory)
    shading_correction(aif_directory,sdc_directory)
    registration_batch(sdc_directory,rgs_directory,raw_directory,ref_cyc=1)
    #registration_separate(sdc_directory,rgs_directory,raw_directory)
    tile_x,tile_y,tiles = get_tiles(raw_directory)
    patch_missing_tile(rgs_directory,tiles)
    template_stitch(rgs_directory,stc_directory,tile_x,tile_y)
    stitch_from_meta(stc_directory,rgs_directory,tile_width=2304) # change to autodetect in the future
    if segment:
        segment_save(stc_directory,dest_directory)

if __name__ == "__main__":
    main('20210928_Combinatorial_Fluorescent_Barcode_2',segment=False)
