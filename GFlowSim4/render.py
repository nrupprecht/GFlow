#!/usr/local/bin/python3

import cv2
import argparse
import os

# Construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-ext", "--extension", required=False, default='bmp', help="extension name. default is 'bmp'.")
ap.add_argument("-o", "--output", required=False, default='output.mp4', help="output video file")
ap.add_argument("-dir", "--directory", required=False, default='RunData/Pos', help="Directory from which we gather images.")
args = vars(ap.parse_args())

# Arguments
ext = args['extension']
output = args['output']
dir_path = args['directory']

print "Looking for images in ", dir_path

images = []
for f in os.listdir(dir_path):
    if f.endswith(ext):
        images.append(f)

if len(images)==0 : 
    print "No images.\n"
    exit()

# Sorted the list of files
# Sorting comes from <https://stackoverflow.com/questions/33159106/sort-filenames-in-directory-in-ascending-order>
images.sort(key=lambda f: int(filter(str.isdigit, f)))

# Determine the width and height from the first image
image_path = os.path.join(dir_path, images[0])
frame = cv2.imread(image_path)
cv2.imshow('video',frame)
height, width, channels = frame.shape

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v') # Be sure to use lower case
out = cv2.VideoWriter(output, fourcc, 20.0, (width, height))

for image in images:

    image_path = os.path.join(dir_path, image)
    frame = cv2.imread(image_path)

    out.write(frame) # Write out frame to video

    cv2.imshow('video',frame)
    if (cv2.waitKey(1) & 0xFF) == ord('q'): # Hit `q` to exit
        break

# Release everything if job is finished
out.release()
cv2.destroyAllWindows()

print("The output video is {}".format(output))
