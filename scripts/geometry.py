"""Centroid Voronoi and Delaunay geometry using the original scientific algorithms."""
import numpy as np
from scipy.spatial import Voronoi, Delaunay, cKDTree

def polygon_area(poly: np.ndarray) -> float:
    if poly.shape[0] < 3:
        return 0.0
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * float(np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))

def clip_polygon_to_rect(poly: np.ndarray, xmin: float, xmax: float, ymin: float, ymax: float) -> np.ndarray:

    def clip_edge(points, inside_fn, intersect_fn):
        if len(points) == 0:
            return points
        out = []
        prev = points[-1]
        prev_in = inside_fn(prev)
        for curr in points:
            curr_in = inside_fn(curr)
            if curr_in:
                if not prev_in:
                    out.append(intersect_fn(prev, curr))
                out.append(curr)
            elif prev_in:
                out.append(intersect_fn(prev, curr))
            (prev, prev_in) = (curr, curr_in)
        return np.array(out, dtype=float)
    pts = np.array(poly, dtype=float)
    pts = clip_edge(pts, lambda p: p[0] >= xmin, lambda p1, p2: np.array([xmin, p1[1] + (p2[1] - p1[1]) * (xmin - p1[0]) / (p2[0] - p1[0] + 1e-15)]))
    pts = clip_edge(pts, lambda p: p[0] <= xmax, lambda p1, p2: np.array([xmax, p1[1] + (p2[1] - p1[1]) * (xmax - p1[0]) / (p2[0] - p1[0] + 1e-15)]))
    pts = clip_edge(pts, lambda p: p[1] >= ymin, lambda p1, p2: np.array([p1[0] + (p2[0] - p1[0]) * (ymin - p1[1]) / (p2[1] - p1[1] + 1e-15), ymin]))
    pts = clip_edge(pts, lambda p: p[1] <= ymax, lambda p1, p2: np.array([p1[0] + (p2[0] - p1[0]) * (ymax - p1[1]) / (p2[1] - p1[1] + 1e-15), ymax]))
    return pts

def voronoi_finite_polygons_2d(vor: Voronoi, radius: float=1000000.0):
    if vor.points.shape[1] != 2:
        raise ValueError('Voronoi points must have two coordinates')
    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    all_ridges = {}
    for ((p1, p2), (v1, v2)) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))
    for (p1, region_idx) in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if all((v >= 0 for v in region)):
            new_regions.append(region)
            continue
        ridges = all_ridges[p1]
        new_region = [v for v in region if v >= 0]
        for (p2, v1, v2) in ridges:
            if v1 >= 0 and v2 >= 0:
                continue
            v = v1 if v1 >= 0 else v2
            tangent = vor.points[p2] - vor.points[p1]
            tangent /= np.linalg.norm(tangent) + 1e-15
            normal = np.array([-tangent[1], tangent[0]])
            midpoint = (vor.points[p1] + vor.points[p2]) / 2
            direction = np.sign(np.dot(midpoint - center, normal)) * normal
            far_point = vor.vertices[v] + direction * radius
            new_vertices.append(far_point.tolist())
            new_region.append(len(new_vertices) - 1)
        vs = np.array([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        ang = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = [v for (_, v) in sorted(zip(ang, new_region))]
        new_regions.append(new_region)
    return (new_regions, np.array(new_vertices, dtype=float))

def voronoi_areas(points_um: np.ndarray, bbox_um: tuple[float, float, float, float]) -> np.ndarray:
    (xmin, xmax, ymin, ymax) = bbox_um
    n = len(points_um)
    if n < 4:
        return np.full(n, np.nan)
    span = max(xmax - xmin, ymax - ymin)
    radius = 10 * span if span > 0 else 1000000.0
    vor = Voronoi(points_um)
    (regions, vertices) = voronoi_finite_polygons_2d(vor, radius=radius)
    areas = np.empty(n, dtype=float)
    for (i, reg) in enumerate(regions):
        poly = vertices[reg]
        poly = clip_polygon_to_rect(poly, xmin, xmax, ymin, ymax)
        areas[i] = polygon_area(poly)
    return areas

def delaunay_neighbor_distance(points_um: np.ndarray) -> float:
    if len(points_um) < 4:
        return np.nan
    tri = Delaunay(points_um)
    edges = set()
    for simp in tri.simplices:
        (i, j, k) = simp
        edges.add(tuple(sorted((i, j))))
        edges.add(tuple(sorted((j, k))))
        edges.add(tuple(sorted((i, k))))
    edges = np.array(list(edges), dtype=int)
    vec = points_um[edges[:, 0]] - points_um[edges[:, 1]]
    d = np.linalg.norm(vec, axis=1)
    return float(np.mean(d)) if len(d) else np.nan

def knn_mean_distance(points_um, k=3):
    """Mean of each centroid's distances to its k nearest other centroids."""
    if len(points_um) < k + 2:
        return np.nan
    distances, _ = cKDTree(points_um).query(points_um, k=k + 1)
    return float(distances[:, 1:k + 1].mean(axis=1).mean())
