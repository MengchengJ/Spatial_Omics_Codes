import os
import re

def run_ijm(file_dir):
    os.system(r'C:\"Program Files"\fiji-win64\Fiji.app\ImageJ-win64.exe' + f' --headless -macro --mem=32000M "{file_dir}"')

def create_mist_command(tile_x,tile_y,input_dir,output_dir,glob_file=None):
    with open('MIST_test.ijm') as f:
        t = f.read()
        t = re.sub('imagedir=.* filenamepattern=',f'imagedir={input_dir} filenamepattern=',t)
        t = re.sub('outputpath=.* display',f'outputpath={output_dir} display',t)
        t = re.sub('gridwidth=.* gridheight=.* starttile',f'gridwidth={tile_x} gridheight={tile_y} starttile',t)
        t = re.sub('extentwidth=.* extentheight=.* timeslices=',f'extentwidth={tile_x} extentheight={tile_y} timeslices=',t)
    with open('MIST_temp.ijm','w') as f:
        f.write(t)
    
if __name__ == "__main__":
    #create_mist_command(5,6,'something_in','something_out')
    run_ijm('MIST_temp.ijm')