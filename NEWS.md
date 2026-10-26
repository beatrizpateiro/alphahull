# alphahull 3.0.0

- Replaced the `interp` triangulation with `geometry::delaunayn` (Qhull) in
  `delvor()`. Delaunay edge adjacencies and Voronoi edges now come from the
  same triangulation, with a convex-hull consistency check.
- Oriented unbounded Voronoi rays using the sample centroid and extended their
  finite plotting markers beyond the midpoint projection when necessary.
- Updated `plot.delvor()` to draw the Delaunay mesh directly from `mesh`.
- **Compatibility change:** `delvor(...)$tri.obj` is now an integer matrix with
  three columns of point indices, rather than an `interp` `triSht` object.
  
## Development notes

* Parts of the code and documentation were revised with the assistance of AI tools.