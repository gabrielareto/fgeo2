# Example datasets

This folder contains real forest census and environmental data, derived quadrat
locations, and artificial examples for surveying and coordinate recovery. The
inventory below includes the files in `coordinate_recovery_test/`. Row counts
refer to the copies stored here and exclude text-file headers.

## Barro Colorado Island censuses, taxonomy and environment

The environmental extraction tutorial uses the BCI 50-ha forest dynamics plot,
with local coordinates spanning 1,000 m west to east and 500 m south to north.

| File | Contents |
| --- | --- |
| `census1.txt` | Existing CTFS tab-separated stem records for census 1, nominally 1982. There are 263,934 rows, with measurement dates from March 1981 to July 1983. |
| `census2.txt` | Existing CTFS tab-separated stem records for census 2 in 1985. There are 299,447 rows. |
| `bci.spptable.rdata` | Historical CTFS species reference: object `bci.spptable`, a data frame with 1,053 rows and 10 columns. |
| `bci-elevation.csv` | Complete historical BCI elevation grid with 20,301 rows and columns `x`, `y`, `z`. Coordinates and elevation are in metres. |
| `bci-soils-wolf-2015.csv` | Complete Wolf et al. soil dataset with 625 sampling locations and all 30 original columns. |

### Tree censuses

The two census files are the historical CTFS text exports already present in
this repository. Their original release version is undocumented here. The
curated BCI archive describes the study, census history, coordinate conventions
and data dictionaries. Credit the BCI investigators Richard Condit, Rolando
Pérez, Salomón Aguilar, Suzanne Lao, Robin Foster and Stephen P. Hubbell.
[BCI data archive, Condit et al. 2019](https://doi.org/10.15146/5xcp-0d46).

In these exports, `PX` and `PY` are plot coordinates in metres, `Mnemonic` is
the species code, and `StemID` identifies a stem within a census. `TreeID` can
repeat because a tree can have several stems. Preserve `Tag`, `StemTag` and
`Quadrat` as character columns when reading the files, so leading zeros survive.
The pair `{Census, StemID}` uniquely identifies each supplied census record.

### Species reference

`bci.spptable.rdata` loads the data frame `bci.spptable`, with species codes,
Latin names, genera, families, taxonomic authorities and related fields. Join
the census `Mnemonic` to its `sp` column. The stored R object was verified as
identical to the one distributed in ForestGEO's
[AGBfluxes data file](https://github.com/forestgeo/AGBfluxes/blob/master/data/bci.spptable.rda),
although the serialized files differ. See the package's
[table description](https://forestgeo.github.io/AGBfluxes/reference/bci.spptable.html).

This is a regional reference table, so its row count is not the species richness
of the BCI 50-ha plot. The
[Panama census archive](https://doi.org/10.15146/mdpr-pm59) explains the regional
scope of the table named `bci.spptable`. Credit the CTFS/ForestGEO taxonomists and
BCI census team. The original download and taxonomy revision of this historical
copy are not recorded.

### Elevation

`bci-elevation.csv` was converted from `BCI elev.TXT` in Gabriel Arellano's
existing BCI analysis files, with a source-file modification date of 19 April
2017. The conversion changed the delimiter to comma and renamed `elev` to `z`.
All numerical values and rows were retained. The grid contains every combination
of `x = 0, 5, ..., 1000` and `y = 0, 5, ..., 500`.

The values were verified against the published `BCIelev.tsv` in the BCI archive
by Condit et al. Normalizing the historical table to tab-separated numeric text
reproduced the archive file's full SHA-256 checksum,
`0869ea08d2488b1e5c67d6e78201c6d0e938123be74ae190210d411aefd683b4`.
The archive defines `x` and `y` as distances from the west and south plot borders,
and elevation relative to sea level. Its methods refer to Harms et al. (2001),
*Habitat associations of trees and shrubs in a 50-ha neotropical forest plot*.
The 5 m spacing describes the supplied grid; its survey and preprocessing history
determine the underlying accuracy. The Dryad release is **CC0 1.0**.
[BCI data archive](https://doi.org/10.15146/5xcp-0d46) ·
[Published elevation file metadata](https://datadryad.org/api/v2/files/148947).

### Soil chemistry

Wolf, J. A., Hubbell, S. P., Fricker, G. A. and Turner, B. L. (2015),
*Geospatial observations on tropical forest surface soil chemistry*,
Ecology 96, 2313. [Paper DOI](https://doi.org/10.1890/15-0558.1).

This file is an unchanged copy of the deposited `soils.csv`, with a more
descriptive filename. It contains 2010 surface-soil measurements, local plot
coordinates `gx`, `gy`, UTM coordinates, nutrients, carbon, nitrogen and pH.
Keep the extraction methods in the variable names. For example, `M3Ca` is
Mehlich-3 calcium and `Ca` is BaCl2-extractable calcium. The header contains
a leading space in ` utm_n`; the tutorial removes header whitespace on import.
[Variable definitions and quality-control notes](https://esapubs.org/archive/ecol/E096/204/metadata.php).

The publisher's dataset is released under **CC0 1.0**. Its original metadata
also explicitly states that the data have no copyright restrictions. Retain
the source citation and identify any derived or interpolated products.
[Dataset and licence](https://doi.org/10.6084/m9.figshare.3562863) ·
[CC0 terms](https://creativecommons.org/publicdomain/zero/1.0/).

The unchanged file has MD5 `d4deb84fa3a2300097efc07a7a42700c`. Its original
carriage-return line endings are retained. [Original CSV download](https://ndownloader.figshare.com/files/5635584).

## Luquillo censuses and quadrat locations

These data represent the 16-ha Luquillo Forest Dynamics Plot in Puerto Rico,
with a 320 m east-west by 500 m north-south footprint. Credit the Luquillo Forest
Dynamics Plot census team, Luquillo LTER and CTFS/ForestGEO. The official
[LFDP description](https://luquillo.lter.network/forest-dynamics/), under
"LFDP censusing", documents the plot, census sequence and 400 quadrats.

| File | Contents |
| --- | --- |
| `luquillo.full1.rdata` | Object `luquillo.full1`: 125,446 rows and 19 columns, for nominal census 1 (1990). Nonmissing measurement dates span 27 June 1990 to 29 September 1993. |
| `luquillo.full2.rdata` | Object `luquillo.full2`: 125,446 rows and 19 columns, for nominal census 2 (1995). Nonmissing measurement dates span 26 October 1994 to 14 May 1998. |
| `luquillo.quadrats.20x20.csv` | Derived lookup table for all 400 quadrats, with columns `quadrat_id`, `xmin`, `ymin`, `xmax`, `ymax`, in metres. IDs run from `0101` to `1625` in the plot's column-row convention. |

The R files are historical CTFS **full tree tables** already included in fgeo2.
Each has one row per tree. Status codes distinguish living trees (`"A"`), dead
trees (`"D"`) and entries before recruitment (`"P"`). Both files have the same
ordered `treeID` vector. Their total row counts therefore include trees without
a measurement in the corresponding census. For multistem trees, the table
contains the diameter of one stem; it is not an inventory of every stem.
Richard Condit's
[CTFS R formatted data tables](https://repository.si.edu/bitstreams/5c31a4aa-6103-4c6b-8ce6-0d6c716d2a62/download),
pp. 1-2, describes these conventions. The original export versions are not
recorded for these copies.

Load each R file with `load()` to recover the object named above. Coordinates
are stored as `gx`, `gy`. The quadrat table was generated and validated within
fgeo2, and added in
[479935a](https://github.com/gabrielareto/fgeo2/commit/479935a978bd711f175e4cf9953f846bf4fb87fd).
The [quadrat tutorial](../03_spatial/tutorials/quadrats.qmd) reproduces and checks
the table against the census coordinates and naming convention. Read
`quadrat_id` as character to preserve leading zeros.

## Other Panama census plots

`P14.csv` and `P19.csv` are historical CTFS measurement exports from two
one-hectare plots in the Smithsonian's Panama plot network. Their `PlotName`
and numeric `PlotID` values match the
[Smithsonian plot register](https://ctfs.si.edu/Public/Datasets/PanamaTreePlots/PanamaPlotInfo.php):
P14 is plot 112 in Panama province, and P19 is plot 117 in Colon province.

| File | Contents |
| --- | --- |
| `P14.csv` | 5,222 measurement records and 32 columns, covering six censuses (`PlotCensusNumber` 1-6), with measurement dates from 8 November 1996 to 14 September 2015. |
| `P19.csv` | 8,339 measurement records and 32 columns, covering seven censuses (`PlotCensusNumber` 1-7), with measurement dates from 27 December 1996 to 4 April 2014. |

Both files are **tab-separated despite their `.csv` extension**; use
`read.delim()`. They contain plot, tree, stem and measurement identifiers,
species, local coordinates `PX`, `PY`, diameter, measurement height, dates and
status. `DBHID` is unique within each file, whereas `{StemID, CensusID}` is not.
Preserve identifiers and tags as character strings when importing.

Credit Richard Condit, Rolando Pérez, Suzanne Lao, Salomón Aguilar and the
Smithsonian/CTFS field teams. The
[Panama plot-network description](https://ctfs.si.edu/Public/Datasets/PanamaTreePlots/)
provides project attribution and points to Condit et al. (2013) for methods.
The exact download URLs and export versions for these historical copies are
not recorded; the matching site identifiers establish which plots they describe.

## Simulated topography survey

| File | Contents |
| --- | --- |
| `topography-survey.csv` | 336 simulated pairwise elevation differences connecting 181 distinct locations in a 240 x 240 m area. |

This example was developed for fgeo2. It mainly follows a planned 20 m grid,
with additional measurements around unusable points. Columns are `from_point`,
`to_point`, `from_x_m`, `from_y_m`, `to_x_m`, `to_y_m` and
`height_difference_m`. Coordinates and differences are in metres; a positive
difference means the destination is higher than the origin. The
[topography-solver tutorial](../06_topography_habitat/tutorials/topography-solver.qmd)
identifies the survey as simulated and documents its use. The same input is used
in the [terrain tutorial](../06_topography_habitat/tutorials/derive-topographical-variables.qmd).

## Practice field map and clicked coordinates

| File | Contents |
| --- | --- |
| `example-map.jpeg` | Scanned practice map created by Gabriel Arellano: a 10 x 10 m quadrat divided into four 5 x 5 m subquadrats, with control points and numbered example trees. |
| `example-map-clicked-coor.csv` | 35 clicked image locations, with columns `id`, `x`, `y`: nine controls (`A`-`D`, `C0`-`C4`) and 26 numbered example points (`0011`-`0036`). |

These are teaching materials for the
[coordinate-recovery tutorial](../06_topography_habitat/tutorials/coordinate-recovery.qmd).
The JPEG is explicitly labelled as a practice map. The CSV stores image
coordinates, not metres. The tutorial supplies the control points' known plot
coordinates, spanning x = 20-30 m and y = 50-60 m, and demonstrates the conversion.
Read `id` as character so labels retain their leading zeros.

## Coordinate-recovery evaluation files

All files in this section are in `coordinate_recovery_test/`. Gabriel Arellano
created an artificial reference grid, printed it, scanned or photographed the
paper under several conditions, and clicked the resulting images to compare
coordinate-recovery methods. The
[ImageJ and geometry evaluation](../evaluation/5-imagej-and-geometry.qmd),
under "Practical comparison of algorithms", contains the reference-grid
generator, treatment descriptions and analysis.

| File | Contents |
| --- | --- |
| `coordinate_recovery_test_grid.pdf` | Generated reference diagram with 64 known points on an 8 x 8 grid from 48 to 62 m on both axes. Controls A-D delimit a 10 x 10 m quadrat at x = 50-60 m, y = 50-60 m. |

Each row below pairs an image with its clicked coordinates.

| Image file | Coordinates file | Paper treatment and capture |
| --- | --- | --- |
| `test A.pdf` | `points test A.csv` | Untreated paper, scanned. |
| `test B.pdf` | `points test B.csv` | Crumpled and flattened paper, scanned. |
| `test C.pdf` | `points test C.csv` | Crumpled, torn into pieces and reassembled paper, scanned. |
| `test D.pdf` | `points test D.csv` | Untreated paper, photographed by hand with a phone. |
| `test E.pdf` | `points test E.csv` | Crumpled and flattened paper, photographed by hand with a phone. |
| `test F.pdf` | `points test F.csv` | Crumpled, torn into pieces and reassembled paper, photographed by hand with a phone. |

Each CSV has 68 rows and columns `id`, `x`, `y`: four controls
`corner_A`-`corner_D`, followed by `p1`-`p64`. Coordinates are in each image's own
pixel frame, with a bottom-left origin. They were obtained with R's `locator()`;
the clicking procedure is retained in
[coordinate-locator.R](../evaluation/coordinate-locator.R). These artificial
test points and physically distorted paper copies are separate from the
practice field map above.

## Attribution and reuse

Dataset-specific terms apply independently of the licence for the R code.
The CC0 statements above refer to the identified elevation and soil releases.
The original data-specific licences for the older census exports and species
reference have not been verified here; consult the providers' terms for their
intended use. Cite the source datasets, retain the distinction between
observations and generated examples, and identify any subsequent transformations.
