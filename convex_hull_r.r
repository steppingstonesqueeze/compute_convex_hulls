# Convex Hull Computation in R
# Using Graham Scan Algorithm

# Function to compute cross product of vectors OA and OB
cross_product <- function(O, A, B) {
  return((A[1] - O[1]) * (B[2] - O[2]) - (A[2] - O[2]) * (B[1] - O[1]))
}

# Function to compute convex hull using Graham Scan
convex_hull <- function(points) {
  n <- nrow(points)
  if (n < 3) {
    stop("Need at least 3 points to form a convex hull")
  }
  
  # Find the bottom-most point (and left-most in case of tie)
  start_idx <- which.min(points[, 2] + points[, 1] * 1e-10)
  start_point <- points[start_idx, ]
  
  # Sort points by polar angle with respect to start point
  other_points <- points[-start_idx, ]
  angles <- atan2(other_points[, 2] - start_point[2], 
                  other_points[, 1] - start_point[1])
  
  # Sort by angle, then by distance for collinear points
  distances <- sqrt((other_points[, 1] - start_point[1])^2 + 
                   (other_points[, 2] - start_point[2])^2)
  sorted_indices <- order(angles, distances)
  sorted_points <- rbind(start_point, other_points[sorted_indices, ])
  
  # Graham scan
  hull <- list()
  for (i in 1:nrow(sorted_points)) {
    # Remove points that make right turn
    while (length(hull) > 1 && 
           cross_product(hull[[length(hull) - 1]], 
                        hull[[length(hull)]], 
                        sorted_points[i, ]) <= 0) {
      hull <- hull[-length(hull)]
    }
    hull[[length(hull) + 1]] <- sorted_points[i, ]
  }
  
  # Convert to matrix
  hull_matrix <- do.call(rbind, hull)
  
  # Ensure clockwise orientation
  hull_matrix <- ensure_clockwise(hull_matrix)
  
  return(hull_matrix)
}

# Function to ensure clockwise orientation
ensure_clockwise <- function(hull) {
  n <- nrow(hull)
  if (n < 3) return(hull)
  
  # Calculate signed area (positive for counterclockwise, negative for clockwise)
  signed_area <- 0
  for (i in 1:n) {
    j <- if (i == n) 1 else i + 1
    signed_area <- signed_area + (hull[j, 1] - hull[i, 1]) * (hull[j, 2] + hull[i, 2])
  }
  
  # If counterclockwise (positive area), reverse the order
  if (signed_area > 0) {
    hull <- hull[nrow(hull):1, ]
  }
  
  return(hull)
}

# Function to get line segments from hull vertices
get_hull_segments <- function(hull) {
  n <- nrow(hull)
  segments <- data.frame(
    x1 = numeric(n),
    y1 = numeric(n),
    x2 = numeric(n),
    y2 = numeric(n)
  )
  
  for (i in 1:n) {
    j <- if (i == n) 1 else i + 1
    segments[i, ] <- c(hull[i, 1], hull[i, 2], hull[j, 1], hull[j, 2])
  }
  
  return(segments)
}

# Function to plot convex hull
plot_convex_hull <- function(points, hull, title = "Convex Hull") {
  plot(points, pch = 16, col = "blue", main = title,
       xlab = "X", ylab = "Y", asp = 1)
  
  # Draw hull
  hull_closed <- rbind(hull, hull[1, ])
  lines(hull_closed, col = "red", lwd = 2)
  
  # Highlight hull vertices
  points(hull, pch = 16, col = "red", cex = 1.2)
  
  # Add vertex labels
  text(hull[, 1], hull[, 2], labels = 1:nrow(hull), 
       pos = 3, col = "red", font = 2)
}

# Example 1: Simple square
cat("Example 1: Simple square\n")
points1 <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0.5, 0.5), ncol = 2, byrow = TRUE)
colnames(points1) <- c("x", "y")
print("Input points:")
print(points1)

hull1 <- convex_hull(points1)
print("Convex hull vertices (clockwise):")
print(hull1)

segments1 <- get_hull_segments(hull1)
print("Hull line segments:")
print(segments1)

# Example 2: Random points
cat("\nExample 2: Random points\n")
set.seed(42)
points2 <- matrix(runif(20, 0, 10), ncol = 2)
colnames(points2) <- c("x", "y")
print("Input points (first 5):")
print(head(points2, 5))

hull2 <- convex_hull(points2)
print("Convex hull vertices (clockwise):")
print(hull2)

segments2 <- get_hull_segments(hull2)
print("Hull line segments:")
print(segments2)

# Create plots
par(mfrow = c(1, 2))
plot_convex_hull(points1, hull1, "Example 1: Square with Interior Point")
plot_convex_hull(points2, hull2, "Example 2: Random Points")

# Example 3: Collinear points
cat("\nExample 3: Points with collinear vertices\n")
points3 <- matrix(c(0, 0, 2, 0, 4, 0, 4, 2, 2, 2, 0, 2, 1, 0, 3, 0), 
                  ncol = 2, byrow = TRUE)
colnames(points3) <- c("x", "y")
print("Input points:")
print(points3)

hull3 <- convex_hull(points3)
print("Convex hull vertices (clockwise):")
print(hull3)

segments3 <- get_hull_segments(hull3)
print("Hull line segments:")
print(segments3)