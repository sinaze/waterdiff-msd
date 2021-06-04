locmsd: main.c tools.c analyze.c
	gcc -o locmsd main.c tools.c analyze.c -lm -lxdrfile -lblas -L/Users/sina_ze/xdrfile -I/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers -I.
