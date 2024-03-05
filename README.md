# lefer

A small C++ library for drawing evenly-spaced and non-overlapping curves in a flow field (also called of "vector field" in some contexts), using the Jobard and Lefer (1997) algorithm.
This algorithm is thoroughly described in a scientific paper ([Jobard and Lefer 1997](#references)), but you might find
[this article useful too](https://pedro-faria.netlify.app/posts/2024/2024-02-19-flow-even/en/index.html).

![](./images/even_curves2.png)


# How to build it?

This project is built by CMake. You can build the project by running:

```bash
cmake .
make
```

# Calculating curves

The functions from this library calculates all coordinates from each curve you want
to draw. They make sure that the coordinates from each curves does not collide (or overlap)
with the coordinates from other curves.

Very briefly, the idea behind the algorithm, is to draw a curve by walking through
the vector field, and constantly check if we are getting to close from neighbouring
curves. If we do get too close, then, we stop drawing the current curve, and
start to draw a different curve in a different position of the flow field.


# The main API

The core part of the Jobard and Lefer algorithm can be splitted in two parts:

- Drawing non-overlapping curves;
- Drawing evenly-spaced, and also, non-overlapping curves;


This library offers a single function for each part (`lefer::even_spaced_curves()` and `lefer::non_overlapping_curves()`).
So, if you want to draw curves that do not overlap each other, but you do not care about
how much far they are from each other, you probably want to use the `lefer::non_overlapping_curves()` function.
Otherwise, you use the `lefer::even_spaced_curves()`.

Both functions return a `std::vector` of `Curve` objects. Each `Curve` object represents a curve that
was drawn into the flow field.

## References

Jobard, Bruno, and Wilfrid Lefer. 1997. “Creating Evenly-Spaced Streamlines of Arbitrary Density.” In Visualization
    in Scientific Computing ’97, edited by Wilfrid Lefer and Michel Grave, 43–55. Vienna: Springer Vienna.
