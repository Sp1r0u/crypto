CC       = gcc
CFLAGS   = -Wall -L.
LDFLAGS  = -lpbc -lgmp -lm
OBJFILES = main.o EC.o field.o poly.o circuit.o wire.o gate.o CRS.o proof.o

TARGET   = zkSNARK

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
