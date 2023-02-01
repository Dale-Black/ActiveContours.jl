# ActiveContours

Image segmentation via active contours. 

## Goals
- [x] Reimplement basic PDE-based (level set) segmentation algorithm(s) from scikit-image
  - [x] `src/chanvese.jl` ([Orinal Chan-Vese Link](https://github.com/scikit-image/scikit-image/blob/v0.19.2/skimage/segmentation/_chan_vese.py#L175-L347))
- [ ] Modify these algorithms to use SciML's ecosystem where available (ModellingToolkit.jl, MOL.jl, etc)
- [ ] Investigate how to add neural networks to these algorithms
