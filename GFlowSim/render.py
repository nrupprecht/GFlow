#!/usr/local/bin/python3

# Original code taken from <http://tsaith.github.io/combine-images-into-a-video-with-python-3-and-opencv-3.html>

import cv2
import argparse
import os
import sys

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-ext", "--extension", required=False, default='bmp', help="extension name. default is 'bmp'.")
ap.add_argument("-o", "--output", required=False, default='output.mp4', help="output video file")
ap.add_argument("-dir", "--directory", required=False, default='RunData', help="Directory in which the run data is stored.")
args = vars(ap.parse_args())

# Arguments
ext = args['extension']
output = args['output']
dir_path = args['directory']

# Look for images in the 'Pos' subdirectory
load_path = dir_path + "/general/Pos"
save_path = dir_path + "/" + output

print ("Looking for images in {}.".format(load_path))
print ("Will save file in {}.".format(save_path))

images = []
for f in os.listdir(load_path):
    if f.endswith(ext):
        images.append(f)

if len(images)==0 : 
    print ("No images.")
    exit()


# Sorted the list of files
# Sorting method comes from <https://stackoverflow.com/questions/33159106/sort-filenames-in-directory-in-ascending-order>
#images.sort(key=lambda f: int(filter(str.isdigit, f)))
# I had to modify the lambda after I changed to python3 via anaconda.
images.sort(key=lambda f: int(''.join(list(filter(str.isdigit, f)))) )

# Determine the width and height from the first image
image_path = os.path.join(load_path, images[0])
frame = cv2.imread(image_path)
#cv2.imshow('video',frame)
height, width, channels = frame.shape

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
out = cv2.VideoWriter(save_path, fourcc, 20.0, (width, height)) ##

# Variables for the progress bar
count = 1
barlength = 25

# Print the start of the progress bar
sys.stdout.write("["+" "*barlength+"]")
sys.stdout.flush()

# Go through all the images
for image in images:

    image_path = os.path.join(load_path, image)
    frame = cv2.imread(image_path)

    out.write(frame) # Write out frame to video

    #cv2.imshow('video',frame)
    if (cv2.waitKey(1) & 0xFF) == ord('q'): # Hit `q` to exit
        break
    
    length = int(25*count/len(images))
    sys.stdout.write("\r[" + "-" * length + " "*(barlength-length) + "]")
    sys.stdout.flush()
    count += 1

print ("")

# Release everything if job is finished
out.release()
cv2.destroyAllWindows()

print ("The output video is {}".format(save_path)) ##
