# linAlg

A collection of numerical matrix decomposition methods, as discussed on [thoughts-on-coding.com](https://thoughts-on-coding.com).

[Numerical Methods in C++ Part 4: Decomposition of Linear Equation Systems](https://thoughts-on-coding.com/2019/06/12/numerical-methods-with-cpp-part-4-introduction-into-decomposition-methods-of-linear-equation-systems)
- LU-Decomposition
- LU-Decomposition with relative scaled pivot strategy
- Cholesky-Decomposition

## Getting Started

To get it up and running you just need to execute:
- `~\linAlg\build\cmake .. -DCMAKE_TOOLCHAIN_FILE={YOUR_PATH_TO_VCPKG}/scripts/buildsystems/vcpkg.cmake`
- `~\linAlg\build\cmake --build . --config Release`

You can execute the program by `./linAlg`

### Prerequisites/Dependencies

- [fmt](http://fmtlib.net/latest/index.html) External library used for formatting and printing results
- [doctest](https://github.com/onqtam/doctest) Feature-rich C++11/14/17/20 single-header testing framework for unit tests and TDD

## Authors

* **Benjamin Mahr** - [Ben1980](https://github.com/Ben1980)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
