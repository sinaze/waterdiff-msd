locmsd: main.c tools.c analyze.c
	gcc -o locmsd main.c tools.c analyze.c -lm -lxdrfile -L/Users/sina_ze/xdrfile -I.
