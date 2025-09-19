import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Union

class ConvexHull:
    """
    A class to compute convex hull of 2D points using Graham Scan algorithm.
    Returns vertices in clockwise order with line segments.
    """
    
    @staticmethod
    def cross_product(O: np.ndarray, A: np.ndarray, B: np.ndarray) -> float:
        """Compute cross product of vectors OA and OB."""
        return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0])
    
    @staticmethod
    def polar_angle_distance(point: np.ndarray, reference: np.ndarray) -> Tuple[float, float]:
        """Calculate polar angle and distance from reference point."""
        dx, dy = point[0] - reference[0], point[1] - reference[1]
        angle = np.arctan2(dy, dx)
        distance = np.sqrt(dx*dx + dy*dy)
        return angle, distance
    
    @classmethod
    def compute_hull(cls, points: Union[List[List[float]], np.ndarray]) -> np.ndarray:
        """
        Compute convex hull using Graham Scan algorithm.
        
        Args:
            points: Array of 2D points [[x1, y1], [x2, y2], ...]
            
        Returns:
            Array of hull vertices in clockwise order
        """
        points = np.array(points)
        n = len(points)
        
        if n < 3:
            raise ValueError("Need at least 3 points to form a convex hull")
        
        # Find the bottom-most point (leftmost if tie)
        start_idx = np.lexsort((points[:, 0], points[:, 1]))[0]
        start_point = points[start_idx]
        
        # Get other points and sort by polar angle
        other_points = np.delete(points, start_idx, axis=0)
        
        # Calculate polar angles and distances
        polar_data = [cls.polar_angle_distance(point, start_point) for point in other_points]
        
        # Sort by angle, then by distance for collinear points
        sorted_indices = sorted(range(len(polar_data)), 
                              key=lambda i: (polar_data[i][0], polar_data[i][1]))
        
        sorted_points = np.vstack([start_point, other_points[sorted_indices]])
        
        # Graham scan
        hull = []
        for point in sorted_points:
            # Remove points that make a right turn
            while (len(hull) > 1 and 
                   cls.cross_product(hull[-2], hull[-1], point) <= 0):
                hull.pop()
            hull.append(point)
        
        hull = np.array(hull)
        
        # Ensure clockwise orientation
        return cls._ensure_clockwise(hull)
    
    @staticmethod
    def _ensure_clockwise(hull: np.ndarray) -> np.ndarray:
        """Ensure hull vertices are in clockwise order."""
        n = len(hull)
        if n < 3:
            return hull
        
        # Calculate signed area using shoelace formula
        signed_area = 0
        for i in range(n):
            j = (i + 1) % n
            signed_area += (hull[j, 0] - hull[i, 0]) * (hull[j, 1] + hull[i, 1])
        
        # If counterclockwise (positive area), reverse the order
        if signed_area > 0:
            hull = hull[::-1]
        
        return hull
    
    @staticmethod
    def get_hull_segments(hull: np.ndarray) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
        """
        Get line segments from hull vertices.
        
        Returns:
            List of line segments as ((x1, y1), (x2, y2)) tuples
        """
        n = len(hull)
        segments = []
        
        for i in range(n):
            j = (i + 1) % n
            segments.append(((hull[i, 0], hull[i, 1]), (hull[j, 0], hull[j, 1])))
        
        return segments
    
    @staticmethod
    def plot_hull(points: np.ndarray, hull: np.ndarray, title: str = "Convex Hull"):
        """Plot the points and their convex hull."""
        plt.figure(figsize=(8, 6))
        
        # Plot all points
        plt.scatter(points[:, 0], points[:, 1], c='blue', s=50, alpha=0.7, label='Points')
        
        # Plot hull vertices
        plt.scatter(hull[:, 0], hull[:, 1], c='red', s=100, label='Hull vertices', zorder=5)
        
        # Draw hull edges
        hull_closed = np.vstack([hull, hull[0]])  # Close the hull
        plt.plot(hull_closed[:, 0], hull_closed[:, 1], 'r-', linewidth=2, label='Hull edges')
        
        # Add vertex numbers
        for i, (x, y) in enumerate(hull):
            plt.annotate(f'{i+1}', (x, y), xytext=(5, 5), textcoords='offset points',
                        fontweight='bold', color='red')
        
        plt.title(title)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        plt.show()


def main():
    """Demonstrate convex hull computation with examples."""
    
    # Example 1: Simple square with interior point
    print("Example 1: Simple square with interior point")
    points1 = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0.5, 0.5]])
    print(f"Input points:\n{points1}")
    
    hull1 = ConvexHull.compute_hull(points1)
    print(f"Convex hull vertices (clockwise):\n{hull1}")
    
    segments1 = ConvexHull.get_hull_segments(hull1)
    print("Hull line segments:")
    for i, segment in enumerate(segments1):
        print(f"  Segment {i+1}: {segment[0]} -> {segment[1]}")
    
    ConvexHull.plot_hull(points1, hull1, "Example 1: Square with Interior Point")
    
    # Example 2: Random points
    print("\nExample 2: Random points")
    np.random.seed(42)
    points2 = np.random.uniform(0, 10, (12, 2))
    print(f"Input points (first 5):\n{points2[:5]}")
    
    hull2 = ConvexHull.compute_hull(points2)
    print(f"Convex hull vertices (clockwise):\n{hull2}")
    
    segments2 = ConvexHull.get_hull_segments(hull2)
    print(f"Number of hull segments: {len(segments2)}")
    print("Hull line segments:")
    for i, segment in enumerate(segments2):
        print(f"  Segment {i+1}: ({segment[0][0]:.2f}, {segment[0][1]:.2f}) -> "
              f"({segment[1][0]:.2f}, {segment[1][1]:.2f})")
    
    ConvexHull.plot_hull(points2, hull2, "Example 2: Random Points")
    
    # Example 3: Collinear points
    print("\nExample 3: Points with collinear vertices")
    points3 = np.array([[0, 0], [2, 0], [4, 0], [4, 2], [2, 2], [0, 2], [1, 0], [3, 0]])
    print(f"Input points:\n{points3}")
    
    hull3 = ConvexHull.compute_hull(points3)
    print(f"Convex hull vertices (clockwise):\n{hull3}")
    
    segments3 = ConvexHull.get_hull_segments(hull3)
    print("Hull line segments:")
    for i, segment in enumerate(segments3):
        print(f"  Segment {i+1}: {segment[0]} -> {segment[1]}")
    
    ConvexHull.plot_hull(points3, hull3, "Example 3: Rectangle with Collinear Points")
    
    # Example 4: Triangle
    print("\nExample 4: Simple triangle")
    points4 = np.array([[0, 0], [4, 0], [2, 3], [1, 1], [2, 1]])
    print(f"Input points:\n{points4}")
    
    hull4 = ConvexHull.compute_hull(points4)
    print(f"Convex hull vertices (clockwise):\n{hull4}")
    
    segments4 = ConvexHull.get_hull_segments(hull4)
    print("Hull line segments:")
    for i, segment in enumerate(segments4):
        print(f"  Segment {i+1}: {segment[0]} -> {segment[1]}")
    
    ConvexHull.plot_hull(points4, hull4, "Example 4: Triangle with Interior Points")


if __name__ == "__main__":
    main()
