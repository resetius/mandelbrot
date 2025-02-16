#include <vector>
#include <cstdint>
#include <string>
#include <cstring>
#include <cmath>
#include <chrono>
#include <map>
#include <random>

#include "big_float.h"

template<typename T>
class MandelbrotRenderer {
private:
    int width, height;
    int max_iterations;
    int frame = 0;
    T center_x, center_y, view_width;
    T eps = 1.0;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;

    void write_iterations(const std::vector<int>& buffer, const std::string& filename) {
        FILE* fp = fopen(filename.c_str(), "wb");
        if (!fp) return;
        fwrite(&width, sizeof(width), 1, fp);
        fwrite(&height, sizeof(height), 1, fp);
        fwrite(buffer.data(), sizeof(int), buffer.size(), fp);
        int blocks = 0;
        if constexpr(!std::is_same_v<T,double>) {
            T tmp;
            blocks = tmp.getMantissa().size();
        }
        fwrite(&blocks, sizeof(blocks), 1, fp);
        fwrite(&center_x, sizeof(center_x), 1, fp);
        fwrite(&center_y, sizeof(center_y), 1, fp);
        fwrite(&view_width, sizeof(view_width), 1, fp);
        fclose(fp);
    }

    std::pair<T, T> find_zoom_point(const std::vector<int>& iterations_buffer) {
        // Находим точку с максимальным градиентом в центральной области
        int center_x = width / 2;
        int center_y = height / 2;
        int search_radius = width / 4;
        double _1_w = 1.0 / width;
        double _1_h = 1.0 / height;

        std::multimap<double,std::pair<int,int>> top_scores;
        int scores = 10;

        // Преобразование экранных координат в координаты комплексной плоскости
        T pixel_width = view_width * T(_1_w);
        T pixel_height = (view_width * T(height) * T(_1_w)) * T(_1_h);

        for (int y = center_y - search_radius; y < center_y + search_radius; y++) {
            for (int x = center_x - search_radius; x < center_x + search_radius; x++) {
                if (x < 2 || x >= width-2 || y < 2 || y >= height-2) continue;

                // Проверяем окрестность точки
                int black_count = 0;
                int colored_count = 0;
                int total_count = 0;

                for (int dy = -2; dy <= 2; dy++) {
                    for (int dx = -2; dx <= 2; dx++) {
                        int idx = ((y + dy) * width + (x + dx));
                        if (iterations_buffer[idx] == max_iterations) {
                            black_count++;
                        } else {
                            colored_count++;
                        }
                        total_count ++;
                    }
                }

                // Ищем области на границе множества
                // double score = black_count * (total_count-colored_count);
                //double score = colored_count == 0
                //    ? 0
                //    : (double)black_count / (double)colored_count;
                double score = black_count == 0
                    ? 0
                    : (double)colored_count / (double)black_count;

                // Штраф за удаление от текущего центра
                double distance_penalty = std::sqrt(
                    std::pow(x - center_x, 2) +
                    std::pow(y - center_y, 2)
                ) / search_radius;
                score *= (1.0 - 0.3 * distance_penalty);

                top_scores.insert(std::make_pair(score, std::make_pair(x, y)));

                if ((int)top_scores.size() > scores) {
                    top_scores.erase(top_scores.begin());
                }
            }
        }

        double r = dist(gen);  // равномерное число в [0,1)
        double exponent = 2.0;
        int n = (int)top_scores.size();

        // Вычисляем индекс с bias:
        // При r = 0 получим index = n-1 (конец), при r = 1 — index = 0 (начало)
        size_t index = static_cast<size_t>( std::floor( (n - 1) - std::pow(r, exponent) * (n - 1) ) );

        // Итерация к нужному элементу (multimap не поддерживает random access)
        auto it = top_scores.begin();
        std::advance(it, index);
        auto [best_x, best_y] = it->second;

        // Преобразуем экранные координаты в координаты множества Мандельброта
        T new_x = this->center_x + T(best_x - width/2) * pixel_width;
        T new_y = this->center_y + T(best_y - height/2) * pixel_height;

        return {new_x, new_y};
    }

public:
    MandelbrotRenderer(
        int w, int h, int max_iter,
        T center_x = -0.75, T center_y = 0, T view_width = 4.0)
        : width(w), height(h), max_iterations(max_iter),
          center_x(center_x), center_y(center_y), view_width(view_width),
          gen(rd()), dist(0.0, 1.0)
    {
        T one(1.0);
        while (one + eps > one) {
            eps = T(0.5) * eps;
        }
        eps = T(4096) * eps;
    }

    std::tuple<T,T,T,int,int> get_parameters() {
        return {center_x,center_y,view_width,frame,max_iterations};
    }

    void render_animation(int from, int num_frames) {
        std::vector<int> iterations_buffer(width * height);
        double _1_w = 1.0 / width;

        for (frame = from; frame < num_frames; frame++) {
            T pixel_size = view_width * T(_1_w);

            auto start = std::chrono::high_resolution_clock::now();

            std::cerr << "Zoom: " << (double)view_width << ", " << max_iterations << "\n";

            #pragma omp parallel for collapse(2)
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    T x0 = center_x + T(x - width/2) * pixel_size;
                    T y0 = center_y + T(y - height/2) * pixel_size;

                    int iter = get_iteration_mandelbrot(x0, y0);

                    iterations_buffer[y * width + x] = iter;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            double seconds = std::chrono::duration<double>(end - start).count();

            std::cerr << "frame: " << frame << ", elapsed: " << seconds << "\n";

            std::string filename = "frame_" + std::to_string(frame) + ".dat";
            write_iterations(iterations_buffer, filename);

            // Находим следующую точку для зума
            auto [new_x, new_y] = find_zoom_point(iterations_buffer);
            center_x = new_x;
            center_y = new_y;

            // Уменьшаем ширину обзора вдвое
            //view_width = view_width * T(0.5);
            view_width = view_width * T(1.0/1.5);
            //max_iterations += 10;
            max_iterations += 1;

            if (view_width < eps) {
                std::cerr << "stop on eps: " << (double)eps << "\n";
                break;
            }
        }
    }

    int get_iteration_mandelbrot(T x0, T y0) {
        T x = 0;
        T y = 0;
        T xn;
        T x2=x*x, y2=y*y;
        T _4 = T(4);
        int i;
        for (i = 1; i < max_iterations && x2+y2 < _4; i=i+1) {
            x2=x*x; y2=y*y;
            xn = x2 - y2 + x0;
            if constexpr(std::is_same_v<T,double> || std::is_same_v<T,_Float128>) {
                y = T(2.0)*x*y + y0;
            } else {
                y = (x*y).Mul2() + y0;
            }
            x = xn;
        }

        return i;
    }

};

constexpr int MAX_PRECISION = 32;

template <int Precision>
struct RenderType {
    using type = BigFloat<Precision, uint64_t>;
};

template <>
struct RenderType<0> {
    using type = double;
};

template <int Precision>
void render_recursive(int width, int height, int frames,
                      int current_frame, int max_iterations,
                      const typename RenderType<Precision>::type& center_x,
                      const typename RenderType<Precision>::type& center_y,
                      const typename RenderType<Precision>::type& view_width)
{
    using T = typename RenderType<Precision>::type;
    MandelbrotRenderer<T> renderer(width, height, max_iterations, center_x, center_y, view_width);
    renderer.render_animation(current_frame, frames);
    auto [nx, ny, nv, nframe, nits] = renderer.get_parameters();

    if (nframe < frames) {
        if constexpr (Precision < MAX_PRECISION) {
            render_recursive<Precision + 1>(width, height, frames, nframe, nits, nx, ny, nv);
        } else {
            std::cerr << "Maximum Precision reached: "
                << Precision << ", frames left "
                << nframe << " of " << frames << std::endl;
        }
    }
}

void render_continue(const std::string& fn, int frames) {
    // fn must be like `frame_number.dat`
    FILE* fp = fopen(fn.c_str(), "rb");
    if (!fp) {
        std::cerr << "Cannot open file " << fn << "\n";
        return;
    }
    size_t prefix_pos = fn.find("frame_");
    if (prefix_pos == std::string::npos) {
        std::cerr << "Filename " << fn << " does not contain 'frame_'\n";
        fclose(fp);
        return;
    }
    prefix_pos += 6; // len of "frame_"
    size_t dot_pos = fn.find('.', prefix_pos);
    if (dot_pos == std::string::npos) {
        std::cerr << "Filename " << fn << " has no extension separator\n";
        fclose(fp);
        return;
    }
    std::string frame_number_str = fn.substr(prefix_pos, dot_pos - prefix_pos);
    int frameno = 0;
    try {
        frameno = std::stoi(frame_number_str);
    } catch (const std::exception& e) {
        std::cerr << "Cannot parse frame number from filename " << fn << "\n";
        fclose(fp);
        return;
    }

    int width;
    int height;
    fread(&width, sizeof(width), 1, fp);
    fread(&height, sizeof(height), 1, fp);
    std::vector<int> iterations;
    iterations.resize(width*height);
    fread(iterations.data(), sizeof(int), width*height, fp);
    int max_iterations = *std::max_element(iterations.begin(), iterations.end());
    int blocks = 0;
    fread(&blocks, sizeof(blocks), 1, fp);

    std::cerr << "Run from checkpoint: "
        << width << " " << height << " " << max_iterations << " "
        << blocks << "\n";

    if (blocks == 0) {
        double center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading double parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<0>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 1) {
        BigFloat<1, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<1,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<1>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 2) {
        BigFloat<2, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<2,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<2>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 3) {
        BigFloat<3, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<3,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<3>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 4) {
        BigFloat<4, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<4,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<4>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 5) {
        BigFloat<5, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<5,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<5>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 6) {
        BigFloat<6, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<6,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<6>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 7) {
        BigFloat<7, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<7,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<7>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 8) {
        BigFloat<8, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<8,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<8>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 9) {
        BigFloat<9, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<9,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<9>(width, height, frames, frameno,
                            max_iterations,
                            center_x, center_y, view_width);
    } else if (blocks == 10) {
        BigFloat<10, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<10,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<10>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 11) {
        BigFloat<11, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<11,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<11>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 12) {
        BigFloat<12, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<12,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<12>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 13) {
        BigFloat<13, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<13,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<13>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 14) {
        BigFloat<14, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<14,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<14>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 15) {
        BigFloat<15, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<15,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<15>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else if (blocks == 16) {
        BigFloat<16, uint64_t> center_x, center_y, view_width;
        if (fread(&center_x, sizeof(center_x), 1, fp) != 1 ||
            fread(&center_y, sizeof(center_y), 1, fp) != 1 ||
            fread(&view_width, sizeof(view_width), 1, fp) != 1) {
            std::cerr << "Error reading BigFloat<16,uint64_t> parameters\n";
            fclose(fp);
            return;
        }
        render_recursive<16>(width, height, frames, frameno,
                             max_iterations,
                             center_x, center_y, view_width);
    } else {
        std::cerr << "Unsupported blocks value: " << blocks << "\n";
    }

    fclose(fp);
}

int main(int argc, char** argv) {
    int frames = 20; //1000;
    int max_iterations = 50; //100;
    int width = 1920;
    int height = 1080;
    std::string restore_point;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--width") && i < argc-1) {
            width = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--height") && i < argc-1) {
            height = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--frames") && i < argc-1) {
            frames = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--max_iterations") && i < argc-1) {
            max_iterations = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--restore") && i < argc-1) {
            restore_point = argv[++i];
        }
    }

    if (!restore_point.empty()) {
        render_continue(restore_point, frames);
    } else {
        render_recursive<0>(width, height, frames, 0, max_iterations, -0.75, 0, 4.0);
    }

    return 0;
}
