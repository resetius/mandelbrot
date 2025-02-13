All: mandelbrot_bigfloat.exe mandelbrot.exe mandelbulb.exe julia.exe mandelbrot_raw.exe mandelbrot_dat2png.exe

clean:
		rm -f *.o *.exe 

mandelbulb.exe: mandelbulb_gui.o mandelbulb.o
		$(CC) $^ $(CFLAGS) `pkg-config --libs gtk+-3.0 epoxy` -lm -o $@

%.exe: %.o Makefile
		$(CXX) $< $(CFLAGS) `pkg-config --libs gtk+-3.0 epoxy libpng` -lm -o $@

%.o: %.c Makefile
		$(CC) -g -Wall $(CFLAGS) `pkg-config --cflags gtk+-3.0 epoxy` -c $< -o $@

%.o: %.cpp Makefile
		$(CXX) -g -Wall -std=c++20 $(CFLAGS) `pkg-config --cflags gtk+-3.0 epoxy` -c $< -o $@

