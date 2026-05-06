# Variablen für den Compiler und die Flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall

TARGET = matrix_calc

SRCS = matrix.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

# Befehl zum Aufräumen: 'make clean'
clean:
	rm -f $(TARGET)