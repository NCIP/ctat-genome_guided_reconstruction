#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build genomeguidedreco.v${VERSION}.simg docker://trinityctat/genomeguidedreco:${VERSION}
