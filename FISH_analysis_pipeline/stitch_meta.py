import os
import re
import numpy as np
import cv2
import pandas as pd
from tqdm import tqdm
from skimage.io import imread
from skimage.io import imsave

def get_meta_df(stc_directory):
    file_location = os.path.join(stc_directory,[f for f in os.listdir(stc_directory) if 'global' in f][0])
    raw = ['\n'.join(line.strip(';').split('; ')) for line in open(file_location).read().strip('\n').split('\n')]
    cols, vals = zip(*[line.split(': ') for line in '\n'.join(raw).split('\n')])
    series = pd.Series(vals,cols)
    series.index = [series.groupby(level=0).cumcount(), series.index]
    meta_df = series.unstack()
    return meta_df

def stitch_from_meta(stc_directory,rgs_directory,tile_width=2304,strip_length=20):
    meta_df = get_meta_df(stc_directory)
    # print([f for f in os.listdir(stc_directory) if '.tif' in f])
    # sample_image_name = [f for f in os.listdir(stc_directory) if '.tif' in f][0]
    # sample_image = imread(os.path.join(stc_directory,sample_image_name))
    meta_df['match'] = meta_df['position'].apply(lambda x: re.match(r'\((\d+)\, *(\d+)\)',x))
    meta_df['y'] = meta_df['match'].apply(lambda x: int(x.group(2)))
    meta_df['x'] = meta_df['match'].apply(lambda x: int(x.group(1)))
    height = max(meta_df['y']) + tile_width
    width = max(meta_df['x']) + tile_width
    stitched_imgs = os.listdir(stc_directory)
    for cyc_chn in tqdm([f for f in os.listdir(rgs_directory) if f.startswith('cyc_')]):
        if f'{cyc_chn}.tif' in stitched_imgs:
            continue
        canvas = np.zeros((height,width),dtype='uint16') # sample_image.shape
        for _,row in meta_df.iloc[::-1].iterrows():
            temp_image = imread(os.path.join(rgs_directory,cyc_chn,row.file))
            # match = re.match(r'\((\d+)\, *(\d+)\)', row.position)
            # x = int(match.group(1))
            # y = int(match.group(2))
            x = row.x
            y = row.y
            canvas[y+strip_length:y+tile_width-strip_length,x+strip_length:x+tile_width-strip_length] = temp_image[strip_length:tile_width-strip_length,strip_length:tile_width-strip_length]
        imsave(os.path.join(stc_directory,f'{cyc_chn}.tif'),canvas,check_contrast=False)
        # cv2.imwrite(os.path.join(stc_directory,f'{cyc_chn}.tif'),canvas)

if __name__ == "__main__":
    stitch_from_meta(r'D:\FISH_images_processed\20210817-sensitivity-03(10um)-quest_processed\stitched',r'D:\FISH_images_processed\20210817-sensitivity-03(10um)-quest_processed\registered','not_important')
                    
