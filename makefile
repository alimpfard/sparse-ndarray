all:
	$(CC) -shared -fpic -c spndarray.c  
	$(CC) -shared -fpic -c spndgetset.c 
	$(CC) -shared -fpic spndarray.o spndgetset.o -o libspndarray.so 

test: all
	$(CC) test.c -L . -lm -lspndarray -o test

clean:
	rm -rf *.so *.o test
