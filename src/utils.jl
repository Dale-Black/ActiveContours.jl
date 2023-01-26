"""
    normalize(img)
Normalize an image based on the mean and the standard deviation

# TODO: check this
"""
normalize(img) = (img .- minimum(img)) ./ (maximum(img) .- minimum(img))
