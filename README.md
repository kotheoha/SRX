# SRX
Here, are presented some indicative code snippets from the project/paper **SRX: Efficient Management of Spatial RDF Data**.
Notice that SRX is an extended/enriched version of popular RDF-3X system (source: `RDF-3X: A RISC-style engine for RDF`).

Below, I mention the lines related to my code contributions/additions relative to the original RDF-3X code (https://code.google.com/archive/p/rdf3x/downloads). I stress that I do that mention **only** for some basic code files; there is a plenty of RDF-3X code files that are affected/enriched by my SRX-contributions in RDF-3X system, so below I just indicatively present how some basic RDF-3X code files are affected/enriched by mine development for SRX purposes.

---

The main component files for the efficient management of `spatial RDF queries` in SRX are:

`Selection.hpp` (path: `rdf3x-version/include/rts/operator/`)
* lines: 36-217, 972-1012

`Selection.cpp` (path: `rdf3x-version/rts/operator/`)
* lines: 31-71, 3270-3322, 3327-3337, 3348-5892, 5974-6981, 7032-7178, 7290-7298, 7311-7315, 7325-7329, 7339-7345

---

The main component files for the efficient management of `spatial RDF updates` in SRX are:

`BulkOperation.hpp` (path: `rdf3x-version/include/rts/runtime/`)
* lines: 9-10, 31-36, 45-65, 80-82, 87-90, 93-94

`BulkOperation.cpp` (path: `rdf3x-version/rts/runtime/`)
* lines: 103-1117, 1350-1397
