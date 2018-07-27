CC       = g++
CFLAGS   = -Wall -L.
LDFLAGS  = -lgmp -lm -lboost_program_options
OBJFILES = main.o field.o element.o

TARGET   = math

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~

