# Recursive Least Squares in Rust

This library implements the Recursive Least Squares algorithm with exponential forgetting.
See [Haykin's Adaptive Filter Theory](http://www.isbnsearch.org/isbn/9780132671453) for
details.

It is using `ndarray` for its vector and matrix data structures.

# Recent releases

+ 0.1.2: the inverse correlation matrix is now correctly initialized as δ^{-1} · 𝟙
+ 0.1.1: use `ndarray`s `Zip`/`NdProducer` functionality via the `azip!` macro for performance;
+ 0.1.0: initial release.


# License

Dual-licensed under Apache 2.0 and MIT licenses to be compatible with the Rust project.
