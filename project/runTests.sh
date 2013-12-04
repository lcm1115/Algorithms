# File: runTests.sh
# Author: Liam Morris
# Description: Runs LCS tests for lcs.cpp implementation.

MAXMEM_MB=2048
ulimit -v $(($MAXMEM_MB * 1024))
n=1
while [ -f args/args.$n ]; do
        cur=1
	args="`cat args/args.$n`"
	echo "lcs $args"
	./lcs $args > out/out.$n
        gprof lcs gmon.out > out/gmon.$n
	n=`expr $n + 1`
done

for f in out/*; do
    if [ $(stat -c%s "$f") == 0 ]; then
        rm "$f"
    fi
done
