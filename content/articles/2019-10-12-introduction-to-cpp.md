Title: Introduction on Numerical Modelling with Differential Equations
date: 2019-10-23 10:00
Modified: 2019-10-23 22:12
comments: true
slug: cpp-101
tags: cpp, linux

<!-- PELICAN_BEGIN_SUMMARY -->
In this section, we summarize the introduction class of C++. Let's go!
<!-- PELICAN_END_SUMMARY -->

C++ is a language both high level and low level, in hence, could be efficient.

This language take extern language, like import in Python, with help of preprocessor, we use `iostream`, `cmath`, `vector` and so on.

In a C++ program always is specified by the keyword return.

In C++ every program has a main function. By UNIX tradition, the error code zero means the action completed successfully. Some either non-zero value means false. In UNIX systems zero represent a message of exit successful.

`.cc` and `.cpp` are valid extensions. In this school we use `.cc` extension for source code and `.hh` for headers files.


### Pure python version ###


``` cpp
// helloworld.cc
#include <iostream>
#include <string>

void print(std::string msg)
{
	std::cout << msg << std::endl;
}

int main(int argc, char** argv)
{
	std::string greeting = "Hello World";
	print(greeting);

	return 0;
}
```

In order to run this program in the terminal, execute this, its depend of your compiler, we will show with gcc, cling, clang and intel compiler.

``` sh
$ g++ -o helloWorld helloworld.cc
$ ./helloworld
```

``` sh
$ clang++ helloworld.cc -o helloworld  
$ ./helloworld
```

``` sh
$ icc helloworld.cc
$ ./a.out
```

``` sh
$ cling
****************** CLING ******************
* Type C++ code and press enter to run it *
*             Type .q to exit             *
*******************************************
[cling]$ #include <iostream>
[cling]$ std::cout << "Hello World" << std::endl;
Hello World
[cling]$ .q
```

Ahora veamos usando CMake

	In [1]: import numpy as np

	In [2]: from memview_bench_v1 import pairwise

	In [3]: X = np.random.random((500, 3))

	In [4]: timeit pairwise(X)
	1 loops, best of 3: 6.51 s per loop

It takes nearly seven seconds to compute 250,000 distances.  This is much
too slow.
<!-- This is a comment -->
{% notebook downloads/notebooks/ChutesAndLadders.ipynb cells[2:] %}