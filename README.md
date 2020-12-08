# Pangenome evaluation

**Experimental and in development.**

## Docker container

For now, feel free to use `jmonlong/job-vgamb:3`.

It was built using [Dockerfile](Dockerfile).
Now that the Snakemake pipeline uses Singularity, the plan is to strip down most tools from it to leave only Snakemake, Singularity and a few basic tools (e.g. bcftools).
The different tools in called in the Snakefile can then used specific docker container (see how seqwish/smoothxg are called in the current Snakefile).
That will make it much easier to define and update the version of specific tools used by the pipeline.

## How do I add and evaluate a new pangenome?

Add your pangenome in GFA/vg/pg format in `s3://vg-k8s/vgamb/chr20`, in a folder following the naming convention: `{method}/{parameters}/{dataset}.{method}.{parameters}.{ext}`.

- The current `{dataset}` we are using is `hpp60-hg38`.
- For example, for a cactus graph: `cactus/{parameters}/hpp60-hg38.cactus.{parameters}.pg`
- No `.` in `{parameters}`, use `-` or `_` instead!

### Short-read mapping evaluation

Run:

```
snakemake --configfile snakemake_config.yaml --config exp=cactus/minigraph-2020-10-28-dna-brnn-nu1-b100k-q5-caf html_out=evaluation-report-srmap.html --cores 8 eval_srmap
```

## I named my graph differently and it doesn't follow the naming convention defined above...

You can update the `input.graphs.s3.links.csv` file.
The first column contains the names of the available file, the second column what it should be called if it followed the name convention. 
When the snakemake pipeline can't find a pangenome, it will look at this file (and create all subsequent file following the usual naming convention).

TL;DR Add a line to `input.graphs.s3.links.csv` with `CURRENT_PATH,CORRECT_PATH`

## Chromosome

To increase turn around time, we do the analysis on one chromosome. 
We started with `chr20` which is the default.
To use a different chromosome, change the [config file](snakemake_config.yaml) or add `chr=chr2` in the `--config` part of the command line.

## Snakemake pipeline

The snakemake pipeline currently implements both the pangenome evaluation and some pangenome construction (minigraph, seqwish/smoothxg, VCF approach).
Eventually, we should maybe split this into separate evaluation and construction pipelines

## Running on our kubernetes server

Start an instance with the `jmonlong/job-vgamb:3` docker image.

The instance needs some files for the snakemake pipeline: Snakefile, config, list of regions, etc
Now that we have a github repo, we could use it to host these files.

We could run the command to make the evaluation report on a big instance and it would run any necessary step.
In practice, we can also evaluate each graph in their own instance (e.g. in case one bugs/gets stuck).

To compute the evaluation stats for one graph, start an instance with e.g. 8 cores + 200 Gb mem + 500 Gb disk and run:

```sh 
git clone https://github.com/jmonlong/pangenome-evaluation
cd pangenome-evaluation

snakemake --use-singularity --configfile snakemake_config.yaml --config exp=cactus/new-cactus-parameters html_out=temp.html --cores 8 eval_srmap --forcerun eval_srmap
```

Once all the pangenomes are created, update the `dataset`/`exp` in `snakemake_config.yaml` to make sure they list the ones to include in the report.
Then, to gather the stats and make the evaluation report, launch a small instance, e.g with 1 core + 10 Gb mem + 100 Gb disk and run:

```sh
git clone https://github.com/jmonlong/pangenome-evaluation
cd pangenome-evaluation

snakemake --use-singularity --configfile snakemake_config.yaml --cores 1 eval_srmap --forcerun eval_srmap
```

The visualization report can be ran in parallel, for example on an instance with 4 cores + 50 Gb mem + 100 Gb disk:

```sh
git clone https://github.com/jmonlong/pangenome-evaluation
cd pangenome-evaluation

snakemake --use-singularity --configfile snakemake_config.yaml --cores 4 viz --forcerun viz
```

## Experimenting with a subset of pangenomes and avoiding conflicts

One could copy the config file and use different output paths.
For example, if it's about exploring cactus parameters, using a [`snakemake_config_cactus.yaml`](snakemake_config_cactus.yaml) that defines their own HTML outputs:

```yaml
html_out: "evaluation-report-cactus.html"
viz_html_out: "visualization-report-cactus.html"
```

Then just use this config in the snakemake commands `--configfile snakemake_config_cactus.yaml`
