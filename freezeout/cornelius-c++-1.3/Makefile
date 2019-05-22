#compiler
CC = g++

#compiler flags
# -g adds debug info
# -Wall turns on most warnings
CFLAGS = -g -Wall

#build target
TARGET = example_cornelius

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
clean:
