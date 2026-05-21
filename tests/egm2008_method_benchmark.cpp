#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <mutex>
#include <thread>
#include <vector>

#include "SHG.hpp"

namespace {

struct Point {
    double r_m;
    double phi_deg;
    double lambda_deg;
};

constexpr double kPi = 3.14159265358979323846;

double deg2rad(double deg) {
    return deg * (kPi / 180.0);
}

double norm3(const std::array<double, 3>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

struct BenchResult {
    double total_ms = 0.0;
    double checksum = 0.0;
};

struct BenchmarkSummary {
    int degree = 0;
    int iterations = 0;
    int threads = 0;
    size_t points_per_iter = 0;
    BenchResult explicit_serial;
    BenchResult model_serial;
    BenchResult explicit_parallel;
    BenchResult model_parallel;
};

BenchResult run_explicit_serial(const std::vector<Point>& points,
                                int degree,
                                int iterations,
                                const std::string& legacy_bin) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        for (const auto& p : points) {
            const double phi = deg2rad(p.phi_deg);
            const double lambda = deg2rad(p.lambda_deg);
            const auto g = SHG::g_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
            sink += norm3(g);
            sink += SHG::U_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
        }
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    return {std::chrono::duration<double, std::milli>(t1 - t0).count(), sink};
}

BenchResult run_model_serial(const std::vector<Point>& points,
                             int degree,
                             int iterations,
                             const SHG::SHGModel& model) {
    volatile double sink = 0.0;
    const int l_max = std::min(degree, model.l_max);
    const int m_max = std::min(degree, model.m_max);
    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        for (const auto& p : points) {
            const double phi = deg2rad(p.phi_deg);
            const double lambda = deg2rad(p.lambda_deg);
            const auto g = SHG::acceleration(model, p.r_m, phi, lambda, l_max, m_max);
            sink += norm3(g);
            sink += SHG::potential(model, p.r_m, phi, lambda, l_max, m_max);
        }
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    return {std::chrono::duration<double, std::milli>(t1 - t0).count(), sink};
}

BenchResult run_explicit_parallel(const std::vector<Point>& points,
                                  int degree,
                                  int iterations,
                                  int thread_count,
                                  const std::string& legacy_bin) {
    double sink = 0.0;
    std::mutex sink_mutex;
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    const auto t0 = std::chrono::high_resolution_clock::now();

    for (int tid = 0; tid < thread_count; ++tid) {
        threads.emplace_back([&, tid]() {
            double local = 0.0;
            for (int i = 0; i < iterations; ++i) {
                for (size_t k = 0; k < points.size(); ++k) {
                    const size_t idx = (k + static_cast<size_t>(tid + i)) % points.size();
                    const auto& p = points[idx];
                    const double phi = deg2rad(p.phi_deg);
                    const double lambda = deg2rad(p.lambda_deg);
                    const auto g = SHG::g_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
                    local += norm3(g);
                    local += SHG::U_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
                }
            }
            std::lock_guard<std::mutex> lock(sink_mutex);
            sink += local;
        });
    }

    for (auto& t : threads) {
        t.join();
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    return {std::chrono::duration<double, std::milli>(t1 - t0).count(), sink};
}

BenchResult run_model_parallel(const std::vector<Point>& points,
                               int degree,
                               int iterations,
                               int thread_count,
                               const SHG::SHGModel& model) {
    double sink = 0.0;
    std::mutex sink_mutex;
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    const int l_max = std::min(degree, model.l_max);
    const int m_max = std::min(degree, model.m_max);
    const auto t0 = std::chrono::high_resolution_clock::now();

    for (int tid = 0; tid < thread_count; ++tid) {
        threads.emplace_back([&, tid]() {
            double local = 0.0;
            for (int i = 0; i < iterations; ++i) {
                for (size_t k = 0; k < points.size(); ++k) {
                    const size_t idx = (k + static_cast<size_t>(tid + i)) % points.size();
                    const auto& p = points[idx];
                    const double phi = deg2rad(p.phi_deg);
                    const double lambda = deg2rad(p.lambda_deg);
                    const auto g = SHG::acceleration(model, p.r_m, phi, lambda, l_max, m_max);
                    local += norm3(g);
                    local += SHG::potential(model, p.r_m, phi, lambda, l_max, m_max);
                }
            }
            std::lock_guard<std::mutex> lock(sink_mutex);
            sink += local;
        });
    }

    for (auto& t : threads) {
        t.join();
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    return {std::chrono::duration<double, std::milli>(t1 - t0).count(), sink};
}

std::string make_legacy_bin(const SHG::SHGModel& model) {
    const std::filesystem::path tmpdir = std::filesystem::temp_directory_path() / "shg_benchmark";
    std::filesystem::create_directories(tmpdir);
    const std::filesystem::path legacy_bin = tmpdir / "EGM2008_legacy.bin";
    if (!SHG::write_coefficients_binary(legacy_bin.string(), model.C, model.S, model.l_max, model.m_max)) {
        throw std::runtime_error("Failed to create legacy bin at " + legacy_bin.string());
    }
    return legacy_bin.string();
}

void print_csv(std::ostream& out, const BenchmarkSummary& summary) {
    const double serial_calls = static_cast<double>(summary.iterations * static_cast<int>(summary.points_per_iter));
    const double parallel_calls = serial_calls * static_cast<double>(summary.threads);

    out << "mode,method,total_ms,per_call_ms,checksum,degree,iterations,threads,points_per_iter\n";
    out << "serial,explicit," << summary.explicit_serial.total_ms << ','
        << (summary.explicit_serial.total_ms / serial_calls) << ','
        << summary.explicit_serial.checksum << ','
        << summary.degree << ',' << summary.iterations << ',' << summary.threads << ',' << summary.points_per_iter << '\n';
    out << "serial,artifact," << summary.model_serial.total_ms << ','
        << (summary.model_serial.total_ms / serial_calls) << ','
        << summary.model_serial.checksum << ','
        << summary.degree << ',' << summary.iterations << ',' << summary.threads << ',' << summary.points_per_iter << '\n';
    out << "parallel,explicit," << summary.explicit_parallel.total_ms << ','
        << (summary.explicit_parallel.total_ms / parallel_calls) << ','
        << summary.explicit_parallel.checksum << ','
        << summary.degree << ',' << summary.iterations << ',' << summary.threads << ',' << summary.points_per_iter << '\n';
    out << "parallel,artifact," << summary.model_parallel.total_ms << ','
        << (summary.model_parallel.total_ms / parallel_calls) << ','
        << summary.model_parallel.checksum << ','
        << summary.degree << ',' << summary.iterations << ',' << summary.threads << ',' << summary.points_per_iter << '\n';
}

void print_human(std::ostream& out, const BenchmarkSummary& summary) {
    const double serial_calls = static_cast<double>(summary.iterations * static_cast<int>(summary.points_per_iter));
    const double parallel_calls = serial_calls * static_cast<double>(summary.threads);

    out << std::fixed << std::setprecision(3);
    out << "degree=" << summary.degree
        << ", iterations=" << summary.iterations
        << ", threads=" << summary.threads
        << ", points_per_iter=" << summary.points_per_iter << std::endl;
    out << "serial explicit_ms=" << summary.explicit_serial.total_ms
        << ", per_call_ms=" << (summary.explicit_serial.total_ms / serial_calls)
        << ", checksum=" << summary.explicit_serial.checksum << std::endl;
    out << "serial artifact_ms=" << summary.model_serial.total_ms
        << ", per_call_ms=" << (summary.model_serial.total_ms / serial_calls)
        << ", checksum=" << summary.model_serial.checksum << std::endl;
    out << "parallel explicit_ms=" << summary.explicit_parallel.total_ms
        << ", per_call_ms=" << (summary.explicit_parallel.total_ms / parallel_calls)
        << ", checksum=" << summary.explicit_parallel.checksum << std::endl;
    out << "parallel artifact_ms=" << summary.model_parallel.total_ms
        << ", per_call_ms=" << (summary.model_parallel.total_ms / parallel_calls)
        << ", checksum=" << summary.model_parallel.checksum << std::endl;
}

} // namespace

int main(int argc, char* argv[]) {
    int degree = 1080;
    int iterations = 200;
    int thread_count = std::max(4u, std::thread::hardware_concurrency());
    std::string csv_output_path;
    bool csv_mode = false;

    if (argc >= 2) {
        degree = std::stoi(argv[1]);
    }
    if (argc >= 3) {
        iterations = std::stoi(argv[2]);
    }
    if (argc >= 4) {
        thread_count = std::max(1, std::stoi(argv[3]));
    }
    for (int i = 4; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--csv" && i + 1 < argc) {
            csv_mode = true;
            csv_output_path = argv[++i];
        } else if (arg == "--human") {
            csv_mode = false;
        }
    }

    const std::vector<Point> points = {
        {6378137.0, 0.0, 0.0},
        {6378237.0, 38.9072, -77.0369},
        {6378137.0, 45.0, 90.0},
        {6378137.0, -45.0, 0.0},
        {6378137.0, 90.0, 0.0},
        {6388137.0, -33.8688, 151.2093},
        {6378137.0, 35.6762, 139.6503},
        {6478137.0, 10.0, 20.0},
    };

    SHG::set_coefficient_loading_verbose(false);
    SHG::set_workspace_enabled(true);

    SHG::SHGModel model;
    try {
        model = SHG::get_model("EGM2008");
    } catch (const std::exception& e) {
        std::cerr << "ERROR: cannot load EGM2008 artifact model: " << e.what() << std::endl;
        return 1;
    }

    const std::string legacy_bin = make_legacy_bin(model);

    // Warm-up both paths.
    {
        const auto& p = points.front();
        const double phi = deg2rad(p.phi_deg);
        const double lambda = deg2rad(p.lambda_deg);
        (void)SHG::g_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
        (void)SHG::U_EGM2008(p.r_m, phi, lambda, degree, legacy_bin);
        (void)SHG::acceleration(model, p.r_m, phi, lambda, std::min(degree, model.l_max), std::min(degree, model.m_max));
        (void)SHG::potential(model, p.r_m, phi, lambda, std::min(degree, model.l_max), std::min(degree, model.m_max));
    }

    const BenchResult explicit_serial = run_explicit_serial(points, degree, iterations, legacy_bin);
    const BenchResult model_serial = run_model_serial(points, degree, iterations, model);

    SHG::set_workspace_enabled(true);
    const BenchResult explicit_parallel = run_explicit_parallel(points, degree, iterations, thread_count, legacy_bin);
    const BenchResult model_parallel = run_model_parallel(points, degree, iterations, thread_count, model);

    BenchmarkSummary summary;
    summary.degree = degree;
    summary.iterations = iterations;
    summary.threads = thread_count;
    summary.points_per_iter = points.size();
    summary.explicit_serial = explicit_serial;
    summary.model_serial = model_serial;
    summary.explicit_parallel = explicit_parallel;
    summary.model_parallel = model_parallel;

    if (csv_mode) {
        if (!csv_output_path.empty()) {
            std::ofstream csv_file(csv_output_path);
            if (!csv_file.is_open()) {
                std::cerr << "ERROR: cannot open CSV output file: " << csv_output_path << std::endl;
                return 1;
            }
            print_csv(csv_file, summary);
        } else {
            print_csv(std::cout, summary);
        }
    } else {
        print_human(std::cout, summary);
    }

    return 0;
}
