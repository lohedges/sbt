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
of the algorithm on the IPU hardware. A number of events are loaded from
file and processed on individual tiles of the IPU. (Modulo the number
of events, i.e. events are replicated to fill the IPU.) When run, the
program measures timing statistics for data transfer to the IPU, the
algorithm run time, and data transfer from the IPU back to the host.
Throughput mesures (number of events per second) are reported as a raw
value, i.e. the algorithm compute speed only, and as an overall value
that accounts for the _full_ run time, i.e. including data transfer
(which will be optimised in future).

Currently the benchmarks are only applicable to MK2 IPU hardware, since the
size of data buffers have been tuned to the size of the MK2 tile memory.
At present, due to memory restrictions, we are only able to utilise 3 of the
hardware threads on each tile, hence throughput could potentially be doubled if
it is possible to restructure the algorithm, or lower the memory footprint of
the intermediate data structures required by the algorithm.

To run the benchmarks:

```
./search_by_triplet num_ipus
```

where `num_ipus` is in the range 1 to 4. (This defaults to 1.) If no IPU device
is found, an `IPUModel` emulator will be used and the specified number of
IPUs is ignored. (In this situation, timing statistics are meaningless.)
Running the program using 4 IPUs should give output like:

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
Running benchmarks on 4 IPUs...
Results...
  Copy to IPU took: 0.001840 ms
  Algorithm execution took: 0.000029 ms
  Raw events per second: 605445.264478
  Copy from IPU took: 0.000636 ms
  Events per second: 7052.187321
````

For the above run, the 4 MK2 IPUs can process events at a throughput of around
605000 events per second, i.e. 605kHz. Accounting for data transfer to/from
the IPU, the throughput drops to around 7kHz. However, we should be able to
interleave data transfer with compute, e.g. ping-ponging data transfer and
compute between IPUs, or streaming many more events to the IPU exchange to be
processed in larger batches.

[In](In) contrast, the CPU implementation of the Search by Triplet algorithm within
Allen can run at 3.4kHz on an Intel Xeon Silver 4215R (3.20GHz), i.e. the IPU
can process about 2 times as many events per second when data transfer is included.
The beefiest GPU currently used for Allen throughput measurements is the NVIDIA
GeForce RTX 3090, which has an approximate throughput of 230kHz for the _entire_
reconstruction sequence!

The next steps are 1) thinking of ways to reduce the algorithm's memory footprint,
so that we can use _all_ of the hardare threads available on the IPU; 2)
developing strategies to utilise the large amount of IPU exchange memory, i.e.
transferring (potentially in the background) events to the exchange, which
could then be processed in larger batches; and 3) looking at ways to minimise
the data transfer cost by interleaving with compute.

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

