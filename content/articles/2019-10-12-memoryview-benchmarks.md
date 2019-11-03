Title: Clase 1
date: 2019-10-12 10:00
Modified: 2019-10-12 22:12
comments: true
slug: simulating-chutes-and-ladders
tags: simulation, animation

<!-- PELICAN_BEGIN_SUMMARY -->
There was recently a [thread](https://groups.google.com/forum/?fromgroups#!topic/cython-users/8uuxjB_wbBQ[1-25] "cython-users archive")
on cython-users which caught my eye.  It has to do with 
[memoryviews](http://docs.cython.org/src/userguide/memoryviews.html), a new
way of working with memory buffers in cython.

I've been thinking recently about how to do fast
and flexible memory buffer access in cython.  I contributed the
[BallTree](http://scikit-learn.org/stable/modules/generated/sklearn.neighbors.BallTree.html)
implementation for nearest neighbors searching in
[scikit-learn](http://www.scikit-learn.org), and have been actively thinking
about how to make it faster and more flexible, including adding the ability
to specify distance metrics other than euclidean and minkowski.

In order to accomplish this, I'd like to have a set of distance metric
functions which take two vectors and compute a distance.  There would
be many functions with similar call signatures which could then be
plugged into a code that would iterate over a set of vectors and
compute the appropriate distances.

<!-- PELICAN_END_SUMMARY -->

$\sum_{-\infty}p_i(x)=\alef$. Abcd



$$
\def\arraystretch{1.5}
   \begin{array}{c:c:c}
   a & b & c \\ \hline
   d & e & f \\
   \hdashline
   g & h & i
\end{array}
$$

$$
\mathcal{O}(n)
$$

### Pure python version ###

In pure python, the implementation described above might look something
like this:

``` python
# memview_bench_v1.py
import numpy as np

def euclidean_distance(x1, x2):
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    return np.sqrt(np.sum((x1 - x2) ** 2))
```
We could exploit symmetry to reduce the number of computations required, but
we'll skip that step for now: this simple version of the function will give
us a good benchmark for comparison with alternatives below.  Using the
`timeit` magic in ipython, we can learn how fast this implementation is:

    In [1]: import numpy as np

    In [2]: from memview_bench_v1 import pairwise

    In [3]: X = np.random.random((500, 3))

    In [4]: timeit pairwise(X)
    1 loops, best of 3: 6.51 s per loop

It takes nearly seven seconds to compute 250,000 distances.  This is much
too slow.
<!-- -->
{% notebook downloads/notebooks/ChutesAndLadders.ipynb cells[2:] %}