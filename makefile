all:
	$(CC) $(CFLAGS) -shared -fpic -c spndarray.c
	$(CC) $(CFLAGS) -shared -fpic -c spndgetset.c
	$(CC) $(CFLAGS) -shared -fpic -c spndreduce.c
	$(CC) $(CFLAGS) -shared -fpic -c spndop.c
	$(CC) $(CFLAGS) -shared -fpic spndarray.o spndgetset.o spndreduce.o spndop.o -o libspndarray.so 

test: all
	$(CC) $(CFLAGS) test.c -L . -lm -lspndarray -o test

clean:
	rm -rf *.so *.o test
