#!/bin/bash

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
	./target/maph "$j"
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

