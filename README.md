# Search by Triplet: IPU

An implementation of the [Search by Triplet](https://doi.org/10.1109/IPDPSW.2019.00118)
track reconstruction algorithm on the
[Graphcore IPU](https://www.graphcore.ai).

The code in this repository has been adapted from parts of the LHCb
[Allen](https://gitlab.cern.ch/lhcb/Allen) project.

## Compiling

First clone the repository:

```
git clone https://github.com/lohedges/sbt.git
```

Now navigate to the repository directory, source your
[Poplar](https://docs.graphcore.ai/projects/poplar-user-guide/en/latest/index.html)
enable script, and run `make`, e.g.:

```
cd sbt
source /opt/poplar/enable.sh
make
```

After this, you should find the `search_by_triplet` executable in the working
directory.

By default the code is compiled using `g++`. To use a different compiler, pass
the `CXX` argument to `make`, e.g.:

```
make CXX=clang++
```

Depending on which version of the Poplar SDK you are using, you might
need to adjust the CXX11 ABI flags. By default, we compile against the
new version of the ABI. To use the old ABI, e.g. on CentOS 7, simply use
the `ABIFLAG` argument when invoking `make`:

```
make ABIFLAG=0
```

The code requires a `C++17` compliant compiler and is also compiled with `-O3`
optimisation flags by default. To pass additional flags, use the `OPTFLAGS`
argument, e.g.:

```
make OPTFLAGS=-funroll-loops
```

## Tests

The code has been tested against the output of Allen and has been
found to return numerically identical tracks in all cases.

## Benchmarks

The `search_by_triplet` executable can be used to measure the performance
of the algorithm on MK2 IPU hardware. A number of events are loaded from
file and processed on individual tiles of the IPU. (Modulo the number
of events, i.e. events are replicated to fill the IPU.) When run, the
program measures timing statistics for running batches of events on four
IPUs. We choose to test two approaches to processing batches of events:

* We use data streams to send event data from the host the IPUs, process
it, the use additional data streams to send back the results. In this
approach, we can utilise 3 threads per tile on each IPU before exhausting
the available memory.

* To run in larger batches we make use of _remote buffers_ on the IPU, which
reside in the IPU exchange memory. This allows us to store replicates of
buffers up to 256MiB in size, from which we can transfer data to/from
the IPU tiles. Due to this size restriction, we can currently only utilise
a single thread per IPU tile due to the size of the remote buffers needed
for the algorithm.

Using remote buffers we also tested different strategies for arranging the
data transfer and compute. We tried a sequential process where data transfer
and compute are alternated for each of the IPUs simultaneously, and a
ping-pong approach where we alternate data transfer and compute between pairs
of IPUs.

* _Sequential_: All IPUs load data from the exchange, then compute, then
copy back to the exchange. This is repeated `num_batches` times.

* Ping-pong: IPUs 0 and 1 copy data from the exchange while IPUs 2 and 3
compute. IPUs 0 and 1 then compute while IPUs 2 and 3 copy results back
to the exchange. This is repeated `num_batches` times. (Note that the first
and last iterations are different, to account for data not being on the
IPU tiles to begin with.)

After testing the above approaches we have found that, when using remote
buffers, throughput saturates after around 50 batches. In addition, the
sequential method was found to give better performance, with a throughput
of roughly 1.5x that of the ping-pong. For our setup, 50 batches with a
single thread per IPU tile on all IPUs is approximately at the limit of the
available remote buffer memory. Reducing the number of batches to 25 allows
us to use two threads per tile, but this was found to result in slightly
reduced throughput. Regardless of the batch size, it is not possible to use
more than two threads, since this puts us beyond the limit of the available
tile memory.

In addition to throughput measurements, the benchmark program also validates
that the output of the replicates is identical.

To run benchmarks using data streams:

```
./search_by_triplet
```

Running the program should give output like:

```
Scanning event directory...
  data/event01.txt
  data/event02.txt
  data/event03.txt
  data/event04.txt
  data/event05.txt
  data/event06.txt
  data/event07.txt
  data/event08.txt
  data/event09.txt
Reading event files...
  data/event01.txt
  data/event02.txt
  data/event03.txt
  data/event04.txt
  data/event05.txt
  data/event06.txt
  data/event07.txt
  data/event08.txt
  data/event09.txt
Creating graph program...
Using data streams...
Results...
  Copy to IPU took: 0.001840 ms
  Algorithm execution took: 0.000029 ms
  Raw events per second: 603031.029777
  Copy from IPU took: 0.000633 ms
  Events per second: 7058.866838
  Compute time fraction: 1.16 %
Validating output...
Finished!
````

To bencmark using remote buffers:

```
./search_by_triplet true
```

You should see output like:

```
Scanning event directory...
  data/event01.txt
  data/event02.txt
  data/event03.txt
  data/event04.txt
  data/event05.txt
  data/event06.txt
  data/event07.txt
  data/event08.txt
  data/event09.txt
Reading event files...
  data/event01.txt
  data/event02.txt
  data/event03.txt
  data/event04.txt
  data/event05.txt
  data/event06.txt
  data/event07.txt
  data/event08.txt
  data/event09.txt
Creating graph program...
Using remote buffers...
Running benchmarks...
Using remote buffers...
Running benchmarks...
Results...
  Copy to remote buffers took: 0.010830 ms
  Algorithm execution took: 0.002572 ms
  Raw events per second: 114466.929966
  Copy from remote buffers took: 0.003039 ms
  Events per second: 17907.093277
  Compute time fraction: 15.64 %
Validating output...
Finished!
````

For the above run, when using data streams, the 4 MK2 IPUs can process events at
a throughput of around 603000 events per second, i.e. 603kHz. (Note that this
timing does not account for data transfer from the host to IPU exchange and back.)
In contrast, the remote buffer approach, which uses a third of the threads per
tile, has a throughput of around 115kHz. When including the cost of data transfer
to and from the host, the throughput falls to around 7kHz when using data streams,
and around 18kHz when using remote buffers. Processing 25 times the number of
events gives a gain of only 2.5x throughput, implying that the time taken for
compute is still too short relative to the data transfer time.

In contrast, the CPU implementation of the Search by Triplet algorithm within
Allen can run at 3.4kHz on an Intel Xeon Silver 4215R (3.20GHz). The beefiest
GPU currently used for Allen throughput measurements is the NVIDIA GeForce
RTX 3090, which has an approximate throughput of 230kHz for the _entire_
reconstruction sequence, i.e. not just the search by triplet algorithm.

## Data

Some example [event data](data) is provided from the [minbias](https://gitlab.cern.ch/lhcb/Allen/-/tree/master/input/minbias)
set provided with the [Allen](https://gitlab.cern.ch/lhcb/Allen) repository.
These files contain hit data for each of the 26 module pairs for each event.
A module pair record starts with information pertaining to the module pair
itself, i.e. the number of hits, and the z coordinate of each module in the
pair. Following this comes the hit records for the module pair. Each record
contains the x, y, and z position of the hit, along with the phi value in
`int16_t` format, and the module index associated with the hit. Within an
event file, module pair records are separated by newlines.
