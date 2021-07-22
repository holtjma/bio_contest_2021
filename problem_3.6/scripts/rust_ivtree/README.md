## Alternate post-contest implementation (part 2)
I was not very satisfied with my original implementation (basically brute-force) from the competition.
This led me to think about this problem a bit more during spare time, and eventually create a more efficient approach that addresses sub-problems 7-9.
The results are preserved in this Rust crate.

This approach includes:
1. Better algorithmic approach - Uses a interval tree of isoform exons (and separate one for introns) to identify the overlaps of each read with any of the exons and then adds the appropriate scoring to that particular isoform. The final step is just an `argmax(...)` of the scores.
2. Implementation in a faster language (Rust)
3. Local parallelization

To run it (assuming Rust and `cargo` are installed and I/O files are in the same location as before):
```
time cargo run --release --bin rust_ivtree
```

The end result is a program that can run on a laptop and solve the problem in roughly 2 hours.
The benchmark using 7 processes (see code to change this):
```
real	137m54.661s
user	930m14.376s
sys	1m51.763s
```
