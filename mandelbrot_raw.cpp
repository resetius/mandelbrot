#include <vector>
#include <cstdint>
#include <string>
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

                for (int dy = -2; dy <= 2; dy++) {
                    for (int dx = -2; dx <= 2; dx++) {
                        int idx = ((y + dy) * width + (x + dx));
                        if (iterations_buffer[idx] == max_iterations) {
                            black_count++;
                        } else {
                            colored_count++;
                        }
                    }
                }

                // Ищем области на границе множества
                double score = black_count * colored_count;

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
            max_iterations += 2;

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
        int i;
        for (i = 1; i < max_iterations && x*x+y*y < 4; i=i+1) {
            if constexpr(std::is_same_v<T,double> || std::is_same_v<T,_Float128>) {
                xn = x*x - y*y + x0;
                y = T(2.0)*x*y + y0;
                x = xn;
            } else {
                xn = x*x - y*y + x0;
                y = (x*y).Mul2() + y0;
                x = xn;
            }
        }

        return i;
    }

};

int main() {
    //MandelbrotRenderer<double> renderer(800, 600, 100);
    //MandelbrotRenderer<BigFloat<2,uint64_t>> renderer(800, 600, 100);
    //renderer.render_animation(300);
    //MandelbrotRenderer<_Float128> renderer(1920, 1080, 100);
    //MandelbrotRenderer<BigFloat<2,uint64_t>> renderer(1920, 1080, 100);
    //MandelbrotRenderer<BigFloat<1,uint64_t>> renderer(1920, 1080, 100);

    //int frames = 200;
    //int max_iterations = 100;
    //int frames = 10;
    int frames = 1000;
    int max_iterations = 50; //100;
    MandelbrotRenderer<double> renderer1(1920, 1080, max_iterations);
    renderer1.render_animation(0, frames);
    auto [x1, y1, v1, frame1, its1] = renderer1.get_parameters();

    if (frame1 < frames) {
        MandelbrotRenderer<BigFloat<1,uint64_t>> renderer2(1920, 1080, its1, x1, y1, v1);
        renderer2.render_animation(frame1, frames);
        auto [x2, y2, v2, frame2, its2] = renderer2.get_parameters();

        if (frame2 < frames) {
            MandelbrotRenderer<BigFloat<2,uint64_t>> renderer3(1920, 1080, its2, x2, y2, v2);
            renderer3.render_animation(frame2, frames);
            auto [x3, y3, v3, frame3, its3] = renderer3.get_parameters();

            if (frame3 < frames) {
                MandelbrotRenderer<BigFloat<3,uint64_t>> renderer4(1920, 1080, its3, x3, y3, v3);
                renderer4.render_animation(frame3, frames);
                auto [x4, y4, v4, frame4, its4] = renderer4.get_parameters();

                if (frame4 < frames) {
                    MandelbrotRenderer<BigFloat<4,uint64_t>> renderer5(1920, 1080, its4, x4, y4, v4);
                    renderer5.render_animation(frame4, frames);
                    auto [x5, y5, v5, frame5, its5] = renderer5.get_parameters();

                    if (frame5 < frames) {
                        MandelbrotRenderer<BigFloat<5,uint64_t>> renderer6(1920, 1080, its5, x5, y5, v5);
                        renderer6.render_animation(frame5, frames);
                        auto [x6, y6, v6, frame6, its6] = renderer6.get_parameters();

                        if (frame6 < frames) {
                            MandelbrotRenderer<BigFloat<6,uint64_t>> renderer7(1920, 1080, its6, x6, y6, v6);
                            renderer7.render_animation(frame6, frames);
                            auto [x7, y7, v7, frame7, its7] = renderer7.get_parameters();
                        }
                    }
                }
            }
        }
    }

    return 0;
}
