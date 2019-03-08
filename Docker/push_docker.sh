#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/genomeguidedreco:${VERSION}
docker push trinityctat/genomeguidedreco:latest
