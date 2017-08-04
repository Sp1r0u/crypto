CC       = gcc
CFLAGS   = -Wall -L.
LDFLAGS  = -lpbc -lgmp -lm
OBJFILES = toolbox.o element.o field.o main.o
TARGET   = field

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~
