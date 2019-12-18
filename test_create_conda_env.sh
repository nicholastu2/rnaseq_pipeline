#!/usr/bin/env bash

if [[ ! -e $HOME/.conda/envs/brentlab_rnaseq_pipeline ]]; then
   conda create -n brent_pipeline_pipeline
elif [[ -e $HOME/.conda/envs/brentlab_rnaseq_pipeline ]]; then
   source activate brentlab_pipeline_pipline
fi
