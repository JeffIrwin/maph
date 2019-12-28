#!/bin/bash

BTYPE=Release

for arg in "$@" ; do
	#echo $arg
	if [[ "$arg" == "Release" ]] ; then
		BTYPE=Release
	elif [[ "$arg" == "Debug" ]] ; then
		BTYPE=Debug
	else
		echo "Warning: unknown cmd argument '$arg'"
		echo
	fi
done

echo "Using build type $BTYPE"
TARGET=target
mkdir -p $TARGET
pushd $TARGET

cmake .. -DCMAKE_BUILD_TYPE=$BTYPE -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON

cmake --build . --config $BTYPE

# from TARGET
popd

