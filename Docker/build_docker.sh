#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t trinityctat/genomeguidedreco:${VERSION} .
docker build -t trinityctat/genomeguidedreco:latest .
