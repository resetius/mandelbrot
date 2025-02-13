#include <vector>
#include <cstdint>
#include <string>
#include <cmath>
#include <chrono>
#include <iostream>

#include <png.h>

class Writer {
private:
    int width, height;

    void write_png(const std::vector<uint8_t>& buffer, const std::string& filename) {
        FILE* fp = fopen(filename.c_str(), "wb");
        if (!fp) return;

        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) {
            fclose(fp);
            return;
        }

        png_infop info = png_create_info_struct(png);
        if (!info) {
            png_destroy_write_struct(&png, nullptr);
            fclose(fp);
            return;
        }

        if (setjmp(png_jmpbuf(png))) {
            png_destroy_write_struct(&png, &info);
            fclose(fp);
            return;
        }

        png_init_io(png, fp);
        png_set_IHDR(png, info, width, height, 8,
                     PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
        png_write_info(png, info);

        std::vector<png_bytep> row_pointers(height);
        for (int y = 0; y < height; y++) {
            row_pointers[y] = (png_bytep)&buffer[y * width * 3];
        }
        png_write_image(png, row_pointers.data());
        png_write_end(png, nullptr);

        png_destroy_write_struct(&png, &info);
        fclose(fp);
    }

public:
    void run() {
        std::vector<uint8_t> buffer;
        std::vector<int> iterations;
        int frame = 0;
        double prev_min_iterations = 0;
        double prev_max_iterations = 0;
        double alpha = 0.2;
        while (true) {
            std::string filename = "frame_" + std::to_string(frame) + ".dat";
            std::string filenamePng = "frame_" + std::to_string(frame) + ".png";
            FILE* f = fopen(filename.c_str(), "rb");
            std::cerr << filename << "->" << filenamePng << "\n";
            if (!f) {
                break;
            }
            fread(&width, sizeof(width), 1, f);
            fread(&height, sizeof(height), 1, f);
            buffer.resize(width*height*3);
            iterations.resize(width*height);
            fread(iterations.data(), sizeof(int), width*height, f);
            double min_iterations = *std::min_element(iterations.begin(), iterations.end());
            double max_iterations = *std::max_element(iterations.begin(), iterations.end());
            std::cerr << min_iterations << " " << max_iterations << "\n";
            double effective_min_iterations = min_iterations;
            double effective_max_iterations = max_iterations;
            if (prev_min_iterations > 0) {
                effective_min_iterations = alpha * effective_min_iterations + (1.0 - alpha) * prev_min_iterations;
                effective_max_iterations = alpha * effective_max_iterations + (1.0 - alpha) * prev_max_iterations;
            }
            prev_min_iterations = effective_min_iterations;
            prev_max_iterations = effective_max_iterations;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    int iter = iterations[y * width + x];

                    double hue = (double)(iter-effective_min_iterations) / (effective_max_iterations-effective_min_iterations);
                    int idx = (y * width + x) * 3;

                    // HSV to RGB преобразование для более интересной визуализации
                    {
                        double h = 360 * hue;
                        double s = 1.0;
                        double v = iter < max_iterations ? 1.0 : 0.0;

                        double c = v * s;
                        double x = c * (1 - std::abs(std::fmod(h / 60.0, 2) - 1));
                        double m = v - c;

                        double r, g, b;
                        if(h < 60) { r = c; g = x; b = 0; }
                        else if(h < 120) { r = x; g = c; b = 0; }
                        else if(h < 180) { r = 0; g = c; b = x; }
                        else if(h < 240) { r = 0; g = x; b = c; }
                        else if(h < 300) { r = x; g = 0; b = c; }
                        else { r = c; g = 0; b = x; }

                        buffer[idx] = static_cast<uint8_t>((r + m) * 255);
                        buffer[idx + 1] = static_cast<uint8_t>((g + m) * 255);
                        buffer[idx + 2] = static_cast<uint8_t>((b + m) * 255);
                    }
                }
            }
            fclose(f);

            write_png(buffer, filenamePng);
            frame++;
        }
    }
};

int main() {
    Writer w;
    w.run();
}