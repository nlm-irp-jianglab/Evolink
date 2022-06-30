# Evolink

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![PyPI](https://img.shields.io/pypi/v/woltka)
![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/woltka)

**Evolink** is a Phylogenetic Approach for Rapid Identification of Genotype-Phenotype Associations in Large-scale Microbial Data

![Evolink](img/Logo.jpg)


## Contents

- [Overview](#overview)
- [Installation](doc/install.md)
- [Example usage](#example-usage)
- Tutorials
  - [Working with WoL](doc/wol.md), [The OGU analysis](doc/ogu.md), [Sequence alignment](doc/align.md)
- Main workflow
  - [Input files](doc/input.md), [Output files](doc/output.md), [Classification systems](doc/hierarchy.md), [Classification methods](doc/classify.md), ["Coord-match" functional profiling](doc/ordinal.md), [Stratification](doc/stratify.md)
- Table utilities
  - [Collapse](doc/collapse.md), [Coverage](doc/coverage.md), [Normalize](doc/normalize.md), [Filter](doc/filter.md), [Merge](doc/merge
  - [QIIME 2](woltka/q2), [Qiita](doc/qiita.md), [Bowtie2](doc/align.md#alignment-with-bowtie2), [SHOGUN](doc/align.md#the-shogun-protocol), [RefSeq](doc/refseq.md), [GTDB](doc/gtdb.md), [MetaCyc](doc/metacyc.md), [KEGG](doc/kegg.md)
- References
  - [Command-line interface](doc/cli.md), [Test datasets](test), [Computational efficiency](doc/perform.md)
- [Citation](#citation)
- [Contact](#contact)

## Overview

### What does Evolink do

Woltka processes [**alignments**](doc/input.md) -- the mappings of microbiome sequencing data against reference sequences (such as genomes or genes), and [infers the best placement](doc/classify.md) of the queries in a hierarchical [classification system](doc/hierarchy.md). One query could have simultaneous matches in multiple references. Woltka finds the most suitable classification unit(s) to describe the query accordingly the criteria specified by the user. Woltka generates [**profiles**](doc/output.md) (feature tables) -- the abundances of classification units which describe the structure or function of microbial communities.

### What else does Woltka do

Woltka provides several utilities for handling feature tables, including [normalizing](doc/normalize.md) data, [collapsing](doc/collapse.md) a table to higher-level features, calculating feature group [coverage](doc/coverage.md), [filtering](doc/filter.md) features based on per-sample abundance, and [merging](doc/merge.md) tables.

### What does Woltka not do

Woltka does NOT **align** sequences. You need to align your sequencing data (FastQ, etc.) against a reference database (we recommend [WoL](wol.md)) using an aligner of your choice (e.g., [Bowtie2](doc/align.md#alignment-with-bowtie2)). The resulting alignment files can be fed into Woltka.

Woltka does NOT **analyze** profiles. We recommend using [QIIME 2](https://qiime2.org/) for robust downstream analyses of the profiles to decode the relationships among microbial communities and with their environments.

Flowchart of Woltka's main classification workflow:

![Woltka process](doc/img/process.png)
