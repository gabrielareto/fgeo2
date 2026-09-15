############################################################
#### DIGITIZE POINTS FROM A SCANNED PAPER MAP
############################################################

# Author: Gabriel Arellano
# Context: CTFS/fgeo2 ImageJ coordinate workflow
# Purpose: read a scanned map stored as PDF, click all mapped
#          points, and save their image coordinates.
#
# Coordinate origin: bottom-left.

# This uses R locator() instead of imageJ.
# Stop clicking by pressing Esc or using the graphics
# device's Finish action.


############################################################
#### INPUT
############################################################

# personal path -- you update
path = "/Users/gabriela/Dropbox/PERSONAL/LLC/Oikobit\ LLC/projects/code-and-statistics-forestgeo/fgeo2/" 

# test case to run
test = "F" # A, B, C, D, E, F

# build complete paths
subpath = "data_examples/coordinate_recovery_test/"
pdf.path = paste0(path, subpath, "test ", test, ".pdf")
csv.path = paste0(path, subpath, "points test ", test, ".csv")

############################################################
#### LOAD AND PLOT
############################################################

map.image <- magick::image_read(pdf.path, density = "300x300")[1] # [1] for first page
map.raster <- as.raster(map.image)

if(test == "A") YLIM = c(1000, 2700)
if(test == "B") YLIM = c(1000, 2700)
if(test == "C") YLIM = c(1000, 2700)
if(test == "D") YLIM = c(1500, 2500)
if(test == "E") YLIM = c(750, 2500)
if(test == "F") YLIM = c(1000, 2500)

plot(map.raster, interpolate = FALSE, ylim = YLIM)
axis(1)
axis(2)
box()


############################################################
#### CLICK AND SAVE POINTS
############################################################

# Click four corners first, in A, B, C, D order:
# A: bottom left
# B: top left
# C: top right
# D: bottom right

# Then the rest of the points. Because this is a test,
# click them in the same numerical order as plotted.

# Then "Esc" to finish, or the other ways; see ?locator

clicks <- locator(
	n = 100,
	type = "p",
	pch = 19,
	col = "red"
)

# Complement with point IDs, separating explicitly the corners
clicks <- data.frame(id = c(paste0("corner_", LETTERS[1:4]), paste0("p", 1:64)), clicks)
head(clicks, 10)

# Save
write.csv(clicks, csv.path, row.names = FALSE)
