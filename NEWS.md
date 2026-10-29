# alphahull 3.0.0

- Replaced the `interp` triangulation with `geometry::delaunayn` (Qhull) in `delvor()`. Delaunay edge adjacencies and Voronoi edges now come from the same triangulation, with a convex-hull consistency check.
- Oriented unbounded Voronoi rays using the sample centroid and extended their   finite plotting markers beyond the midpoint projection when necessary.
- Updated `plot.delvor()` to draw the Delaunay mesh directly from `mesh`.
- **Backward compatibility:** The structure of `delvor(...)$tri.obj` has been preserved as a `triSht`-compatible object to maintain compatibility with existing code and dependent packages, while Delaunay triangulation is now computed using `geometry`.
  
## Development notes

* Parts of the code and documentation were revised with the assistance of AI tools.