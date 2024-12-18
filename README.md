# parquet2bcf

This project reads variant data from a Parquet file, processes it, and writes the results to a BCF (Binary Call Format) file. 

## Features

- Reads variant data from a Parquet file. There is one row per non-reference genotype.
- Writes the processed data to a BCF file.

## Installation

To build and run this project, you need to have Rust installed. You can install Rust using [rustup](https://rustup.rs/).

Clone the repository and navigate to the project directory, then build the project:

```sh
cargo build --release
```

## Usage

### Command-line Arguments

- `--parquet-path`: Path to the input Parquet file containing variant data (default: `test_variants.parquet`).
- `--samples-path`: Path to the file containing sample names (default: `samples.txt`).
- `--output-path`: Path to the output BCF file (default: `output.bcf`).
- `--n-threads`: Number of threads to use for writing the BCF file (default: `2`).

## Example

```sh
target/release/parquet2bcf --parquet-path variants.parquet \
    --samples-path samples.txt \
    --output-path output.bcf \
    --n-threads 4 # output threads
```

## Logging

The project uses the `log` crate for logging. By default, it logs informational messages. You can adjust the logging level by setting the `RUST_LOG` environment variable:

```sh
RUST_LOG=info cargo run --release
```

## Testing

The project includes unit tests for the BCF writing functionality. To run the tests, use the following command:

```sh
cargo test
```

## Contact
Email Josh Weinstock <josh.weinstock@emory.edu> . 
