#!/usr/bin/env python
# coding: utf-8

# ## Setup ####################################################################
import deeplabcut
import cv2 as cv
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import time
import pandas as pd
from pathlib import Path

# ## Flags ####################################################################
write_tracked_video = False
cleanup_file_names = True

# ## Paths ####################################################################

# path to date directories to analyze
date_folders_base =  'D:\RawData'

###############################################################################
# Run this to check if folders are correct

# Option 1: Get list of all dates in date_folders_base (or supply specific list below)
days_list = []
for file in os.listdir(date_folders_base):
    d = os.path.join(date_folders_base, file)
    if os.path.isdir(d):
        if file.startswith('20'):
            days_list.append(file)
    days_list.sort()

# Option 2: Provide a list of dates
# days_list = ['2021-12-08','2021-12-10', '2021-12-11']

print(days_list)

###############################################################################
###############################################################################
## Edit paths to files
###############################################################################
###############################################################################

# Check for xlrd

# Base path for config files
config_base_path = 'D:\King'
os.chdir(config_base_path)

# deeplabcut fly tracker config.yaml path
config_path = os.path.join(config_base_path, 'fly_tracker_2-king-2021-09-27/config.yaml')

# deeplabcut plate corner finder config.yaml path
plate_config_path = os.path.join(config_base_path, 'plate_finder-kinglab-2021-05-02/config.yaml')

# DLC model string
dlc_model_string = 'DLC_resnet50_fly_tracker_2Sep27shuffle1_1000000'

# Plate finder string
dlc_plate_string = 'DLC_resnet50_plate_finderMay2shuffle1_500000'

###############################################################################

def rmdir(directory):
    directory = Path(directory)
    for item in directory.iterdir():
        if item.is_dir():
            rmdir(item)
        else:
            item.unlink()
    directory.rmdir()


def cleanup_directory(cleanup_dir):
    for root, dirs, files in os.walk(cleanup_dir):
        # splitvideos folders
        print("Removing splitvideos directories")
        for name in dirs:
            if name.endswith('splitvideos'):
                print(os.path.join(root, name))
                rmdir(os.path.join(root, name))

        # labeled, segmentation, h5. pickle files
        print("Removing labeled points movies")
        for name in files:
            if name.endswith('_labeled.mp4'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print("Removing plate segmentation images")
        for name in files:
            if name.endswith('_segmentation.jpg'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print("Removing h5 files")
        for name in files:
            if name.endswith('.h5'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print('\n\n')


def rotate_image(img, angle):
    """
    Rotate an image (np.array) by some arbitrary angle about the center of the image.

    See: https://stackoverflow.com/a/9042907/168137

    Parameters:
    img (np.array): Image as np.array, e.g., from opencv read()
    angle (float): Angle for rotation, negative angles are clockwise.

    Returns:
    np.array: Image of same dimensions as image.

    """

    image_center = tuple(np.array(img.shape[1::-1]) / 2)
    rot_mat = cv.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv.warpAffine(img, rot_mat, img.shape[1::-1], flags=cv.INTER_LINEAR)
    return result


def crop(img, x1, y1, x2, y2):
    """
    Crop an image using the upper-left and lower-right corners of a rectangle.

    Parameters:
    img (np.array): Image as np.array, e.g., from opencv read()
    x1, y1 (int): Upper-left corner
    x2, y2 (int): Lower-right corner

    Returns:
    np.array: Cropped image.

    """

    img_crop = img[y1:y2, x1:x2, :]
    return img_crop


def order_points(pts):
    # https://www.pyimagesearch.com/2014/08/25/4-point-opencv-getperspective-transform-example/

    # initialzie a list of coordinates that will be ordered
    # such that the first entry in the list is the top-left,
    # the second entry is the top-right, the third is the
    # bottom-right, and the fourth is the bottom-left
    rect = np.zeros((4, 2), dtype = "float32")

    # the top-left point will have the smallest sum, whereas
    # the bottom-right point will have the largest sum
    s = pts.sum(axis = 1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]

    # now, compute the difference between the points, the
    # top-right point will have the smallest difference,
    # whereas the bottom-left will have the largest difference
    diff = np.diff(pts, axis = 1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]

    # return the ordered coordinates
    return rect


def four_point_transform(image, pts):
    # obtain a consistent order of the points and unpack them
    # individually
    rect = order_points(pts)
    (tl, tr, br, bl) = rect

    # compute the width of the new image, which will be the
    # maximum distance between bottom-right and bottom-left
    # x-coordiates or the top-right and top-left x-coordinates
    widthA = np.sqrt(((br[0] - bl[0]) ** 2) + ((br[1] - bl[1]) ** 2))
    widthB = np.sqrt(((tr[0] - tl[0]) ** 2) + ((tr[1] - tl[1]) ** 2))
    maxWidth = max(int(widthA), int(widthB))

    # compute the height of the new image, which will be the
    # maximum distance between the top-right and bottom-right
    # y-coordinates or the top-left and bottom-left y-coordinates
    heightA = np.sqrt(((tr[0] - br[0]) ** 2) + ((tr[1] - br[1]) ** 2))
    heightB = np.sqrt(((tl[0] - bl[0]) ** 2) + ((tl[1] - bl[1]) ** 2))
    maxHeight = max(int(heightA), int(heightB))

    # now that we have the dimensions of the new image, construct
    # the set of destination points to obtain a "birds eye view",
    # (i.e. top-down view) of the image, again specifying points
    # in the top-left, top-right, bottom-right, and bottom-left
    # order
    dst = np.array([
        [0, 0],
        [maxWidth - 1, 0],
        [maxWidth - 1, maxHeight - 1],
        [0, maxHeight - 1]], dtype = "float32")

    # compute the perspective transform matrix and then apply it
    M = cv.getPerspectiveTransform(rect, dst)
    warped = cv.warpPerspective(image, M, (maxWidth, maxHeight))

    # return the warped image
    return warped

###############################################################################

# Join the days_list with the date folder base to make full paths to the
# directories to be analyzed.
days_recording = [os.path.join(date_folders_base, x) for x in days_list]
print(days_recording)

# Load plate coordinates. Convert to integers.
pl = pd.read_excel(os.path.join(config_base_path, 'plate_coordinates.xls')).astype('int')

# Iterate through folders
for jj, date in enumerate(days_recording):
    print("Processing " + date +"\n")

    # Cleanup directory
    print("Cleaning up directory")
    cleanup_directory(date)

    # Find .avi files
    files = []
    for r, d, f in os.walk(date):
        f.sort()
        for file in f:
            if file.endswith(".avi"):
                files.append(file)
    print(files)

    for ii, vidFile in enumerate(files):
        # Setup path for split video files
        splVidPath = os.path.join(date_folders_base, date, vidFile[:-4] + '_splitvideos')
        print('Videos will write into: ' + splVidPath)

        # Make directory for split videos
        try:
            os.mkdir(splVidPath)
        except OSError as error:
            print(error)

        # Find corners of plate
        deeplabcut.analyze_videos(plate_config_path,
                                  os.path.join(date_folders_base, date, vidFile),
                                  videotype='.avi',
                                  save_as_csv=False)

        # Write movie with plate corners
        time.sleep(10)
        deeplabcut.create_labeled_video(plate_config_path,
                                  os.path.join(date_folders_base, date, vidFile),
                                  videotype='.avi')

        # Rename plate corner finder output files
        old_file =  os.path.join(date_folders_base, date, vidFile[:-4] + dlc_plate_string + '.h5')
        os.rename(old_file, old_file.replace(dlc_plate_string, ''))

        old_file = os.path.join(date_folders_base, date, vidFile[:-4] + dlc_plate_string + '_labeled.mp4')
        os.rename(old_file, old_file.replace(dlc_plate_string, ''))

        # Remove pickle created by DLC
        os.remove(os.path.join(date_folders_base, date, vidFile[:-4] + dlc_plate_string + '_meta.pickle'))

        # Create videoWriters
        for vw_id in range(34):
            outfile = os.path.join(splVidPath, vidFile[:-4] + '_' + str(vw_id) + '.mp4')
            globals()['vW{}'.format(vw_id)] = cv.VideoWriter(
                filename=outfile,
                apiPreference=cv.CAP_FFMPEG,
                fourcc=cv.VideoWriter_fourcc('m', 'p', '4', 'v'),
                fps=4,
                frameSize=(pl.w[vw_id], pl.h[vw_id]))


        # Processing video ###############################################
        # Load video
        vid = cv.VideoCapture(os.path.join(date_folders_base, date, vidFile))
        nframes = int(vid.get(cv.CAP_PROP_FRAME_COUNT))

        # Load csv with tracked corners
        corners = pd.read_hdf(os.path.join(date_folders_base, date,
                              vidFile[:-4] + '.h5'))
        corners = corners[dlc_plate_string]

        # Extract coordinates
        UL = corners['topleft']
        UR = corners['topright']
        LL = corners['bottomleft']
        LR = corners['bottomright']

        # Iterate through frames
        print("\nSplitting into separate videos and cropping.\n")
        for ii in range(nframes):

            # Load image
            ret, img = vid.read()

            # Corner points
            corner_pts = np.array([
                [LL.iloc[ii, 0], LL.iloc[ii, 1]],
                [UL.iloc[ii, 0], UL.iloc[ii, 1]],
                [UR.iloc[ii, 0], UR.iloc[ii, 1]],
                [LR.iloc[ii, 0], LR.iloc[ii, 1]]
                ])

            # Use four point transformation
            warped = four_point_transform(img, corner_pts)

            # Resize to 1200 x 1000
            dim = (1200, 1000)

            # resize image
            resized = cv.resize(warped, dim, interpolation = cv.INTER_AREA)

            # Plate segmentation check. Use 100th frame.
            if ii == 99:
                seg_check = resized.copy()
                for idx, row in pl.iterrows():
                    cv.rectangle(seg_check, (row.x1, row.y1), (row.x2, row.y2), (255,0,0), 1)
                    cv.putText(seg_check, str(idx), (row.x1+15, row.y1+15),
                               cv.FONT_HERSHEY_SIMPLEX, 0.5, (255,0,0), 1, cv.LINE_AA)
                ret = cv.imwrite(os.path.join(date_folders_base, date, vidFile[:-4] + '_segmentation.jpg'), seg_check)

            ## Segment video and write frame to that videoWriter
            for idx, row in pl.iterrows():
                # Create copy of img
                img_tmp = resized.copy()

                # Crop image
                img_crop = crop(img_tmp, row.x1, row.y1, row.x2, row.y2)
                h,w = img_crop.shape[:2]

                # Write image
                globals()['vW{}'.format(idx)].write(img_crop)

        # Release videoWriters
        for vw_id in range(34):
            globals()['vW{}'.format(vw_id)].release()

        # Release main video file
        vid.release()

        # Track points on cropped videos
        print('Pausing before tracking points.\n')
        time.sleep(10)
        deeplabcut.analyze_videos(config_path,
                                  [splVidPath],
                                  videotype='.mp4',
                                  save_as_csv=True)

        # Write movie with tracked fly
        if write_tracked_video:
            print('\nPausing before writing videos with tracked points.\n')
            time.sleep(10)
            deeplabcut.create_labeled_video(config_path,
                                            [splVidPath],
                                            videotype='.mp4',
                                            trailpoints=20)

        # rename output files
        if cleanup_file_names:
            for r, d, f in os.walk(splVidPath):
                for file in f:
                    if file != None and dlc_model_string in file:
                        os.rename(os.path.join(r, file),
                                  os.path.join(r, file.replace(dlc_model_string, '')))

        # Pause before starting the next video
        print('\nPausing before processing the next video.\n')
        time.sleep(20)
