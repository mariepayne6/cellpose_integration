# Cell segmentation and marker quantification using CellPose
# Marie Payne
# Last Update: 3.22.25

# ----------- User Notes ---------------- #
# 1. Make sure to install the required packages
# 2. Make sure your image does not have other samples in the background you do not want to quantify
# 3. Set parameters

# pip install cellpose, 'cellpose[gui]', "opencv-python-headless<4.3", skikit-image, matplotlib

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from cellpose import models
from tifffile import imwrite
import time
from skimage.measure import regionprops, label
import csv

# ------------------ Set Params ------------------ #
DIAMETER = int(19) #@param {type:"number"}, choosing diameter = 0 means cellpose will estimate
SEG_CHANNEL = int(2) # If the image has only one channel, leave it as 0
model = "cyto3" #@param ["cyto", "nuclei", "cyto2", "tissuenet", "livecell"]
channels = [0,1,3]
# Define background thresholds for each channel (adjust these values)
  # I perform manual cell segmentation on one image to determine the thresholds & check output in range
  # Add in 5% margin of error for sensitive stains
background_thresholds = {
    0: 3500,  # CTIP
    1: 2900,  # TBR2
    2: 800,    # Hoechst
    3: 1200    # PAX6
}
img_dir = "/path/to/dir/" 
# Optional: Enter image extension here to read only files/images of specified extension (.tif,.jpg..): 
image_format = "tif" 
SPEED = "slow" # for mask visualization
on = False # Set to True to visualize mask
# Adjust based on # channels for output csv
col_names = ["Filename", 
          "Channel 0 Total", "Channel 0 Positive", "Channel 0 Percentage",
          "Channel 1 Total", "Channel 1 Positive", "Channel 1 Percentage",
          "Channel 3 Total", "Channel 3 Positive", "Channel 3 Percentage"]

# ------------------ Functions ------------------ #
# Function takes in a multichannel tiff file and returns a mask from CellPose
# Inputs are multichannel tiff, cell diameter, model, and channel to segment
# Outputs are segmentation masks
def run_cellpose_masks(img, DIAMETER, SEG_CHANNEL, model, display): 
    print("Using model: ",model)
    model = models.Cellpose(model_type = model)

    if img.ndim > 2:  # Multi-channel image
      img = img[:, :, SEG_CHANNEL]
    else:  # Single-channel image
      img = img

    # Run segmentation (adjust thresholds if needed)
    start_time = time.time()
    masks, flows, styles, diams = model.eval(
            img, 
            diameter=DIAMETER, 
            flow_threshold=0.4, 
            cellprob_threshold=0, 
            channels=[0, 0])
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Cellpose segmentation completed in {elapsed_time:.2f} seconds")

    # Display image
    if display:
      print("Displaying image")
      plt.figure(figsize=(8, 8))
      plt.imshow(img, cmap='gray')  # Show original image
      plt.imshow(masks, cmap='jet', alpha=0.5)  # Overlay segmentation mask
      plt.axis('off')
      plt.title("Cellpose Segmentation Mask")
      plt.show()

    return masks

# Function takes in an input directory and returns a list of files and images
def import_files(input_dir): 
  # r=root, d=directories, f = files
  files=[]

  for r, d, f in os.walk(input_dir):
      for fil in f:
        if (image_format):
          if fil.endswith(image_format):
            files.append(os.path.join(r, fil))
        else:
          files.append(os.path.join(r, fil))
      break #only read the root directory; can change this to include levels

  if(len(files)==0):
    raise RuntimeError("No images found in directory: %s." % input_dir)
  else:
    print("Number of images loaded: %d." %(len(files)))

  # Read and Load all the images
  imgs=[] #store all images

  for f in files:
    im=io.imread(f)
    print("Image loaded. Shape:", im.shape)
    n_dim=len(im.shape) #shape of image
    dim=im.shape #dimensions of image
    channel=min(dim) #channel will be dimension with min value usually
    channel_position=dim.index(channel)
    #if no of dim is 3 and channel is first index, swap channel to last index
    if n_dim==3 and channel_position==0: 
      im=im.transpose(1,2,0)
      dim=im.shape
    
    imgs.append(im)

  return files, imgs

# Function to save the mask
# Inputs: Mask, Image filename
def save_mask(mask, filename):
    output_mask_filename = os.path.splitext(os.path.basename(filename))[0] + "_mask.tif"
    output_mask_path = os.path.join(save_dir, output_mask_filename)
    print("Output mask saved as: ", output_mask_path)
    imwrite(output_mask_path, (mask > 0).astype(np.uint8) * 255)

# Function to visualize the mask
def visualize_mask(test_mask, stain_threshold, labeled_mask, marker, speed):
  plt.figure(figsize=(12, 5))
  plt.subplot(1, 2, 1) # Show segmentation channel with masks
  #plt.imshow(marker, cmap='gray')
  plt.imshow(test_mask, cmap='jet', alpha=0.5)
  plt.title("Segmentation (Channel %d)" % c)
  plt.axis('off')
  
  if (speed == "fast"):
    plt.subplot(1, 2, 2)
    plt.imshow(marker, cmap='gray')  # Again show the marker channel
    plt.scatter([prop.centroid[1] for prop in properties if prop.mean_intensity > stain_threshold], 
                  [prop.centroid[0] for prop in properties if prop.mean_intensity > stain_threshold], 
                  color='red', label='Positive cells', alpha=0.3)
    plt.title(f"Positive Cells (Channel {c})")
    plt.axis('off')
  elif (speed == "slow"):
    # SLOWER & BETTER CODE
    positive_mask = np.zeros_like(test_mask)  # Initialize an empty mask
    for prop in properties:
        if prop.mean_intensity > stain_threshold:
            positive_mask[labeled_mask == prop.label] = 1  # Mark positive regions
    # Show stain positivity mask in the second subplot
    plt.subplot(1, 2, 2)
    plt.imshow(marker, cmap='gray')  # Show the original marker channel
    plt.imshow(positive_mask, cmap='Reds', alpha=0.5)  # Overlay positive cells in red
    plt.title(f"Positive Cells (Channel {c})")
    plt.axis('off')  
  else: 
    print("Invalid speed option. Choose 'fast' or 'slow'.")

  plt.tight_layout()
  plt.show()
  return

# ------------------ Start Code ------------------ #

# Enter Directory path containing the images: 
Input_Directory = img_dir
input_dir = os.path.join(Input_Directory, "") #adds separator to the end regardless if path has it or not

save_dir = input_dir+"Masks/"

if not os.path.exists(save_dir):
  os.makedirs(save_dir)
print("Directory for Masks: ", save_dir)

files, imgs = import_files(input_dir)
layer_counts = [] # Initialize empty csv for storing meta data

# Run Cellpose on all images in directory
for i in range(len(imgs)): # Indexing starts at 0
  img = imgs[i]
  
  print("Running Cellpose on image %s" % files[i])
  test_mask = run_cellpose_masks(img, DIAMETER, SEG_CHANNEL, model, display=on)
  save_mask(test_mask, files[i])

  # Measure Stain Counts
  if img.ndim == 2:  # Single-channel image
    print("Image is single channel, skipping counts")
    pass
  elif img.ndim == 3:  # Multi-channel image
    row = [files[i]] # Store filename
    for c in channels:
      marker = img[:, :, c]
      # Label connected components in mask
      labeled_mask = label(test_mask)

      # Measure properties of labeled nuclei on stain channel
      properties = regionprops(labeled_mask, intensity_image=marker)

      # Define stain positivity threshold (adjust as needed)
      stain_threshold = background_thresholds.get(c, 100)  # Default to 100 if not specified

      # Count positive cells
      positive_cells = sum(prop.mean_intensity > stain_threshold for prop in properties)
      print(f"Total cells detected: {len(properties)}")
      print("For channel %d, number of positive cells: %d" % (c, positive_cells))

      # For troubleshooting
      #for prop in properties:
        #print(prop.mean_intensity)

      if on:
        visualize_mask(test_mask, stain_threshold, labeled_mask, marker, speed=SPEED)

      if len(properties) > 0:
        percentage = 100 * positive_cells / len(properties)
      else:
        percentage = 0  # Avoid division by zero
      print(f"Percentage of positive cells: {percentage:.2f}%")
    
      row.extend([len(properties), positive_cells, f"{percentage:.2f}"])
    layer_counts.append(row)

# Save the layer counts to a CSV file
print("Writing to CSV file")
csv_filename = os.path.join(save_dir, "layer_counts.csv")
headers = col_names

# Write to CSV file once at the end
with open(csv_filename, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(headers)
    writer.writerows(layer_counts)

print(f"CSV file '{csv_filename}' saved successfully.")

