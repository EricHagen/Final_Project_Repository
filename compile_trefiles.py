#! usr/bin/env python3

import os
import shutil

#Make chunked & unchunked folders on computer to which to copy files
os.mkdir("/home/hagen/Desktop/chunked_folder")
os.mkdir("/home/hagen/Desktop/unchunked_folder")

#Set names for those folders
chunked_folder = "/home/hagen/Desktop/chunked_folder"
unchunked_folder = "/home/hagen/Desktop/unchunked_folder"

#Set path to folder with "best" .tre files on Beaulieu1
folder_path = "/home/hagen/Desktop/matK_Sim_Files"
os.chdir(folder_path)

#Establish substrings to look for, for both chunked & unchunked best trees:
best_tree = "bestTree"
chunked_substring = "chunked"

#Write a for loop to compile .tre files into a .trees file
for filename in folder_path:
	if best_tree in filename:
		if chunked_substring in filename:
			current_name = folder_path + filename
			new_name = chunked_folder + filename
			shutil.copy(current_name, new_name)
		else:
			current_name = folder_path + filename
			new_name = unchunked_folder + filename
			shutil.copy(current_name, new_name)
