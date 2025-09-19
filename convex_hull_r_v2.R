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

# Enhanced function to plot convex hull with better visualization
plot_convex_hull <- function(points, hull, title = "Convex Hull") {
  # Set up plot with nice margins and styling
  par(mar = c(4, 4, 3, 1), bg = "white")
  
  # Create base plot
  plot(points, pch = 16, col = "lightblue", main = title,
       xlab = "X", ylab = "Y", asp = 1, cex = 1.2,
       panel.first = grid(col = "lightgray", lty = "dotted"))
  
  # Add border to all points
  points(points, pch = 1, col = "blue", cex = 1.2, lwd = 1.5)
  
  # Draw hull polygon with fill
  hull_closed <- rbind(hull, hull[1, ])
  polygon(hull_closed, col = rgb(1, 0, 0, 0.1), border = NA)
  
  # Draw hull edges with thicker lines
  lines(hull_closed, col = "red", lwd = 3)
  
  # Highlight hull vertices with larger points
  points(hull, pch = 16, col = "red", cex = 1.8)
  points(hull, pch = 1, col = "darkred", cex = 1.8, lwd = 2)
  
  # Add vertex labels with better positioning
  text(hull[, 1], hull[, 2], labels = 1:nrow(hull), 
       pos = 3, col = "white", font = 2, cex = 1.2, offset = 0.8)
  text(hull[, 1], hull[, 2], labels = 1:nrow(hull), 
       pos = 3, col = "darkred", font = 2, cex = 1.2, offset = 0.7)
  
  # Add legend
  legend("topright", legend = c("All Points", "Hull Vertices", "Hull Edges"),
         col = c("lightblue", "red", "red"), 
         pch = c(16, 16, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 3),
         bg = "white", box.col = "gray")
}

# Function to create comprehensive visualization with multiple subplots
create_hull_analysis <- function(points, title_prefix = "Convex Hull Analysis") {
  hull <- convex_hull(points)
  segments <- get_hull_segments(hull)
  
  # Create 2x2 subplot layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # Plot 1: Basic visualization
  plot_convex_hull(points, hull, paste(title_prefix, "- Overview"))
  
  # Plot 2: Step-by-step construction visualization
  plot(points, pch = 16, col = "lightblue", 
       main = paste(title_prefix, "- Construction Steps"),
       xlab = "X", ylab = "Y", asp = 1, cex = 1.2)
  grid(col = "lightgray", lty = "dotted")
  
  # Show construction order with arrows
  for (i in 1:(nrow(hull)-1)) {
    arrows(hull[i, 1], hull[i, 2], hull[i+1, 1], hull[i+1, 2],
           col = "red", lwd = 2, length = 0.1)
    text(hull[i, 1], hull[i, 2], labels = i, 
         col = "white", font = 2, cex = 1.2)
    text(hull[i, 1], hull[i, 2], labels = i, 
         col = "darkred", font = 2, cex = 1.1)
  }
  # Close the hull
  arrows(hull[nrow(hull), 1], hull[nrow(hull), 2], hull[1, 1], hull[1, 2],
         col = "red", lwd = 2, length = 0.1)
  text(hull[nrow(hull), 1], hull[nrow(hull), 2], labels = nrow(hull), 
       col = "white", font = 2, cex = 1.2)
  text(hull[nrow(hull), 1], hull[nrow(hull), 2], labels = nrow(hull), 
       col = "darkred", font = 2, cex = 1.1)
  
  # Plot 3: Distance analysis
  hull_center <- colMeans(hull)
  distances <- sqrt(rowSums((points - matrix(rep(hull_center, nrow(points)), 
                                           nrow = nrow(points), byrow = TRUE))^2))
  hull_distances <- sqrt(rowSums((hull - matrix(rep(hull_center, nrow(hull)), 
                                               nrow = nrow(hull), byrow = TRUE))^2))
  
  plot(points, pch = 16, col = heat.colors(length(distances))[rank(distances)], 
       main = paste(title_prefix, "- Distance from Center"),
       xlab = "X", ylab = "Y", asp = 1, cex = 1.5)
  grid(col = "lightgray", lty = "dotted")
  points(hull_center[1], hull_center[2], pch = 4, col = "black", cex = 2, lwd = 3)
  hull_closed <- rbind(hull, hull[1, ])
  lines(hull_closed, col = "red", lwd = 3)
  
  # Plot 4: Area and perimeter info
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  # Calculate hull area using shoelace formula
  hull_area <- 0
  n <- nrow(hull)
  for (i in 1:n) {
    j <- if (i == n) 1 else i + 1
    hull_area <- hull_area + hull[i, 1] * hull[j, 2] - hull[j, 1] * hull[i, 2]
  }
  hull_area <- abs(hull_area) / 2
  
  # Calculate perimeter
  hull_perimeter <- sum(sqrt(rowSums((segments[,3:4] - segments[,1:2])^2)))
  
  # Display statistics
  text(0.5, 0.9, paste(title_prefix, "- Statistics"), 
       cex = 1.5, font = 2, adj = 0.5)
  text(0.1, 0.7, paste("Total Points:", nrow(points)), cex = 1.2, adj = 0)
  text(0.1, 0.6, paste("Hull Vertices:", nrow(hull)), cex = 1.2, adj = 0)
  text(0.1, 0.5, paste("Hull Area:", round(hull_area, 3)), cex = 1.2, adj = 0)
  text(0.1, 0.4, paste("Hull Perimeter:", round(hull_perimeter, 3)), cex = 1.2, adj = 0)
  
  # Show segments
  text(0.1, 0.25, "Line Segments (clockwise):", cex = 1.2, font = 2, adj = 0)
  for (i in 1:min(4, nrow(segments))) {  # Show first 4 segments
    seg_text <- paste("Seg", i, ": (", round(segments[i,1], 2), ",", 
                     round(segments[i,2], 2), ") → (", 
                     round(segments[i,3], 2), ",", round(segments[i,4], 2), ")")
    text(0.1, 0.15 - (i-1)*0.03, seg_text, cex = 0.9, adj = 0, family = "mono")
  }
  if (nrow(segments) > 4) {
    text(0.1, 0.05, paste("... and", nrow(segments) - 4, "more segments"), 
         cex = 0.9, adj = 0, style = "italic")
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  return(list(hull = hull, segments = segments, area = hull_area, perimeter = hull_perimeter))
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

# Create comprehensive visualizations for all examples
par(mfrow = c(2, 2))

# Example 1: Enhanced visualization
result1 <- create_hull_analysis(points1, "Example 1: Square")

cat("EXAMPLE 1 RESULTS:\n")
cat("Hull area:", result1$area, "\n")
cat("Hull perimeter:", result1$perimeter, "\n")

# Example 2: Enhanced visualization  
result2 <- create_hull_analysis(points2, "Example 2: Random Points")

cat("EXAMPLE 2 RESULTS:\n")
cat("Hull area:", result2$area, "\n")
cat("Hull perimeter:", result2$perimeter, "\n")


# Interactive-style single plot function
plot_interactive_style <- function(points, hull, title = "Convex Hull") {
  # Set up high-quality plot
  par(mar = c(4, 4, 3, 1), bg = "white")
  
  # Create base plot
  plot(points, pch = 16, col = "lightblue", main = title,
       xlab = "X", ylab = "Y", asp = 1, cex = 1.5,
       panel.first = grid(col = "lightgray", lty = "dotted"))
  
  # Add point borders
  points(points, pch = 1, col = "blue", cex = 1.5, lwd = 1.5)
  
  # Draw hull with semi-transparent fill
  hull_closed <- rbind(hull, hull[1, ])
  polygon(hull_closed, col = rgb(1, 0, 0, 0.15), border = NA)
  
  # Draw hull edges
  lines(hull_closed, col = "red", lwd = 4)
  
  # Hull vertices with double outline for visibility
  points(hull, pch = 16, col = "red", cex = 2.2)
  points(hull, pch = 1, col = "darkred", cex = 2.2, lwd = 3)
  
  # Vertex labels with white outline for readability
  text(hull[, 1], hull[, 2], labels = 1:nrow(hull), 
       pos = 3, col = "white", font = 2, cex = 1.4, offset = 1)
  text(hull[, 1], hull[, 2], labels = 1:nrow(hull), 
       pos = 3, col = "darkred", font = 2, cex = 1.3, offset = 0.9)
  
  # Enhanced legend
  legend("topright", 
         legend = c("Input Points", "Hull Vertices", "Hull Boundary", "Hull Fill"),
         col = c("lightblue", "red", "red", rgb(1, 0, 0, 0.15)), 
         pch = c(16, 16, NA, 15), 
         lty = c(NA, NA, 1, NA), 
         lwd = c(NA, NA, 4, NA),
         bg = "white", 
         box.col = "gray",
         cex = 0.7, pt.cex = 0.8,
         inset = 0.02, bty = "o")
}

# Demo with enhanced single plots
cat("ENHANCED SINGLE PLOT DEMONSTRATIONS:\n")

#par(mfrow = c(2, 2))
plot_interactive_style(points1, hull1, "Example 1: Square + Interior (Enhanced)")
Sys.sleep(3)
plot_interactive_style(points2, hull2, "Example 2: Random Points (Enhanced)")
Sys.sleep(3)
# Bonus: Create a challenging example
set.seed(123)
challenging_points <- rbind(
  matrix(runif(16, 0, 8), ncol = 2),  # Random scattered
  matrix(c(2, 2, 6, 2, 6, 6, 2, 6), ncol = 2, byrow = TRUE)  # Square corners
)
challenging_hull <- convex_hull(challenging_points)
plot_interactive_style(challenging_points, challenging_hull, "Bonus: Mixed Distribution")
Sys.sleep(3)
#par(mfrow = c(1, 1))
