
# Consider a 1000 x 1000 grid of 1m cells
resoln <- 1
grid_mtx <- matrix(NA, nrow = 1000, ncol = 1000)

# Let's say the density of fish per m2 in this grid is 1/(50*50), which
#   is equivalent to one fish per 50m cell in our main study.
# We fill in the number of fish expected in each 5m grid cell.
density_mtx <- grid_mtx
density_mtx[] <- (resoln*resoln) * (1/(50*50)) # area * density
N <- sum(density_mtx)

# Let's establish some other true parameters.
sigma <- exp(3)
p0 <- 1

# Say we have an ROV survey that does a 100m transect whose center runs from
# (500, 450) to (500, 550). It perfectly covers an area of 1000 cells (100 cells long by 10 cells wide).
# By our current logic, the expected count at this ROV would be:
N * 1000 / prod(dim(grid_mtx))


# What if, instead, we consider the expected count as the sum over all cells of
# (fraction of density of the Gaussian space use of a fish centered 
#             in that cell contained in the transect) x (num. fish in that cell)

corners <- matrix(c(495, 450,
                    505, 450,
                    505, 550,
                    495, 550), nrow = 4, byrow = T)

# Function to calculate the pct of a Gaussian contained in a rectangle
# (courtesy of ChatGPT)
gaussian_fraction_in_rotated_rect <- function(corners, a = 0, b = 0, sigma = 1) {
  # corners: 4x2 matrix of rectangle corners (x, y)
  # a, b: center of 2D Gaussian
  # sigma: std dev of isotropic Gaussian
  
  # Step 1: Shift rectangle so Gaussian is at origin
  shifted_corners <- corners - matrix(c(a, b), nrow = 4, ncol = 2, byrow = TRUE)
  
  # Step 2: Get vector along one side (use corners 1 and 2)
  v <- shifted_corners[2, ] - shifted_corners[1, ]
  theta <- atan2(v[2], v[1])  # rotation angle of the rectangle
  
  # Step 3: Build rotation matrix (clockwise by theta)
  R <- matrix(c(cos(-theta), -sin(-theta),
                sin(-theta),  cos(-theta)), ncol = 2)
  
  # Step 4: Rotate all shifted corners
  rotated_corners <- t(R %*% t(shifted_corners))
  
  # Step 5: Get axis-aligned bounding box in rotated space
  x1 <- min(rotated_corners[, 1])
  x2 <- max(rotated_corners[, 1])
  y1 <- min(rotated_corners[, 2])
  y2 <- max(rotated_corners[, 2])
  
  # Step 6: Use closed-form CDF
  fx <- pnorm(x2 / sigma) - pnorm(x1 / sigma)
  fy <- pnorm(y2 / sigma) - pnorm(y1 / sigma)
  
  return(fx * fy)
}

# Test: a point on the corner of the rectangle with small sigma should have value 0.25
gaussian_fraction_in_rotated_rect(corners, 495, 450, sigma = 0.1)

# Test: a point on the edge of the rectangle with small sigma should have value 0.5
gaussian_fraction_in_rotated_rect(corners, 495, 500, sigma = 0.1)

# Test: as sigma gets bigger the fraction should reduce as there's some chance
#       the fish extends N/S/E of the transect
gaussian_fraction_in_rotated_rect(corners, 495, 500, sigma = 10)
gaussian_fraction_in_rotated_rect(corners, 495, 500, sigma = 20)


# Now, what we need to do is iterate over the whole grid. For each cell, we
# want to calculate the percentage of the buffer of a hypothetical fish
# centered on that cell contained in our transect.
overlap_from_here_mtx <- grid_mtx
for (i in 1:nrow(grid_mtx)) {
  for (j in 1:ncol(grid_mtx)) {
    overlap_from_here_mtx[i, j] <- gaussian_fraction_in_rotated_rect(corners, i - 0.5, j - 0.5, sigma = 40)
  }
}

image(overlap_from_here_mtx)

# Finally, the expected count we're after is the sum of these overlaps times
# the expected number of fish.
sum(overlap_from_here_mtx * density_mtx)

# it's exactly double the previous expected count??

# Note: I could precompute all the transformations/shifts for each transect/cell pair
# and provide this to NIMBLE, although the RAM cost may not be worth it.

