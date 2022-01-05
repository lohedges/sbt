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
program measures the _raw_ IPU performance in terms of the number of events
processed per second. This excludes any overhead from compiling the graph
program, loading it, streaming data to/from the IPU, etc., i.e. it is the
raw throughput once data is on the IPU.

Running using 1000 tiles on a Graphcore Collosus Mk1 gives a throughput of
around 34000 events per second.(Memory limitations of the Mk1 hardware means
that we can't utilise all tiles at present.) In contrast, the CPU
implementation of the Search by Triplet algorithm within Allen
can process around 3400 on an Intel Xeon Silver 4215R (3.20GHz), i.e.
the IPU can process about 10 times as many events per second. The beefiest
GPU currently used for Allen throughput measurements is the NVIDIA GeForce
RTX 3090, which has an approximate throughput of 230 kHz for the _entire_
sequence, i.e. 7 tiles faster than the IPU performs for a single algorithm!

An RTX 3090 has 10496 CUDA cores. Due to memory limitations mentioned above,
the IPU implementation is currently only using 1000 tiles, each with a single
thread, i.e. the GPU can utilise around 10x more threads. In addition, the GPU
has 24 GB GDDR6X memory. Each IPU tile has 256 kB of SRAM, so around 300 MB
for an entire MK1 device, i.e .much less than the GPU. (There is, of course,
the sizeable IPU exchange memory, but you need to configure/compile data
transfer to / from the IPU tiles where compute is needed.)

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

