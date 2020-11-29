#!/usr/bin/env bash

TAG=temp
DATA_DIR="$(pwd)/tmp"

if [[ ! -e $DATA_DIR ]]; then
  mkdir -p $DATA_DIR
fi

docker build -f sensitivity.Dockerfile -t $TAG-sensitivity .
docker run -v $DATA_DIR:/data $TAG-sensitivity 
