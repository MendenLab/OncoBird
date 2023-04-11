#!/bin/bash

which="shiny"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FILESIZE=$(wc -c <"$DIR/images/oncobird_$which.tar.gz")
echo $FILESIZE
if [[ -f "$DIR/images/oncobird_$which.tar.gz" && "$FILESIZE" -gt "100000" ]]; then
	echo "Image for OncoBird already exists as environment/images/oncobird_$which.tar.gz"
else 
	echo "Image does not exist, creating ..."
	docker build --no-cache=true -f environment/docker_$which/dockerfile_$which -t "oncobird_$which" .
	echo "Export image ..."
        mkdir -p $DIR/images
	docker export $(docker create "oncobird_$which") | gzip -c > $DIR/images/oncobird_$which.tar.gz
fi;


