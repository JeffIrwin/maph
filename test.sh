#!/bin/bash -v

# This works on Windows too, but it needs the right environment.  Launch a
# developer command prompt for VS, then from within that, launch git bash
# ("C:\Program Files\Git\git-bash.exe").

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac
echo ${machine}

if [[ $machine == "Linux" || $machine == "Mac" ]]; then
	target=./target
else
	target=./target/Release
fi

git submodule update --init --recursive

ls submodules/lodepng/
ls submodules/
ls

./clean.sh
./build.sh

echo "==============================================================================="
echo ""
echo "Running tests..."
echo ""

nfail=0
ntotal=0

for j in ./data/*.json; do

	#d=$(dirname "$j")
	p="${j%.json}.png"

	#echo "j = $j"
	#echo "p = $p"
	#echo "d = $d"
	#echo ""

	ntotal=$((ntotal + 1))
	rm "$p"
	${target}/maph "$j"

	if [[ "$?" != "0" ]]; then
		nfail=$((nfail + 1))
		echo "test.sh:  error:  cannot run test $j"
	fi

	diff ./data/expected-output/$(basename $p) $p > /dev/null

	if [[ "$?" == "1" ]]; then
		nfail=$((nfail + 1))
		echo "test.sh:  error:  difference in $j"
	fi

done

echo ""
echo "==============================================================================="
echo ""
echo "Total number of tests  = $ntotal"
echo "Number of failed tests = $nfail"
echo "Done!"
echo ""

exit $nfail
