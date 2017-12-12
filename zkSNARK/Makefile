CC       = gcc
CFLAGS   = -Wall -L.
LDFLAGS  = -lpbc -lgmp -lm
###OBJFILES = BN.o field.o wire.o gate.o circuit.o poly.o CRS.o main.o
OBJFILES = BN.o field.o wire.o gate.o circuit.o CRS.o proof.o main.o

TARGET   = zkSNARK

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
