default: all

CC = gcc -O3

SRC_DIR = ./src/

OBJ = farsa.o

# link math package
LINK = -lm

all: farsa

farsa: $(SRC_DIR)$(OBJ) $(SRC_DIR)client.c  
	$(CC) $(SRC_DIR)$(OBJ) $(SRC_DIR)client.c -o farsa $(LINK)

clean:
	rm $(SRC_DIR)*.o
	rm farsa