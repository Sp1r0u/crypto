CC       = gcc
CFLAGS   = -Wall -L.
LDFLAGS  = -lpbc -lgmp -lm
###OBJFILES = circuit.o element.o field.o BN.o main.o
##OBJFILES = circuit.o field.o BN.o main.o
OBJFILES = BN.o field.o circuit.o main.o
TARGET   = zkSNARK

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
