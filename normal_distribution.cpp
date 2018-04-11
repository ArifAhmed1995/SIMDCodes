#include <chrono>
#include <ctime>

#include <umesimd/UMESimd.h>

#include </include/vdt/vdtMath.h>

const int SIMD_SIZE = 16;

using umesimd_dv = typename UME::SIMD::SIMDVec<double, SIMD_SIZE>;


void compute_with_avx2() {}
void compute_with_svml() {}


void compute_with_umesimd(umesimd_dv &x, umesimd_dv &mean,
                          umesimd_dv &stddev, umesimd_dv &result)
{
    umesimd_dv inv = 1/stddev;
    umesimd_dv err = umesimd_dv(-0.5f) * ((x - mean) * (x - mean) * inv * inv);
    result = inv * umesimd_dv(0.39894228f) * err.exp(); 
}

void compute_with_vdt(double * x, double * mean,
                      double * stddev, double * result, int sample_size)
{
    double inv, err;
    for(int i = 0;i < sample_size;i++)
    {
        inv = 1/stddev[i];
        err = -0.5f * ((x[i] - mean[i]) * (x[i] - mean[i]) * inv * inv);
        result[i] = inv * 0.39894228f * vdt::fast_exp(err);
    }
}

void compute_with_naive_scalar(double * x, double * mean, double * stddev, double * result, int sample_size)
{
    double inv, err;
    for(int i = 0;i < sample_size;i++)
    {
        inv = 1/stddev[i];
        err = -0.5f * ((x[i] - mean[i]) * (x[i] - mean[i]) * inv * inv);
        result[i] = inv * 0.39894228f * exp(err);
    }
}

int main()
{
    int sample_size = 10000;

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end; 
    double x[sample_size], mean[sample_size], stddev[sample_size], result[sample_size]; 

    for(int i = 0;i < sample_size;i++)
    {
        x[i] = 10.0f + drand48() * 5.0f;
        mean[i] = 3.0f + drand48() * 3.0f;
        stddev[i] = 2.0f + drand48() * 4.0f;
        result[i] = 0.0;
    }

    //UME::SIMD computation.
    umesimd_dv x_vec(x);
    umesimd_dv mean_vec(mean);
    umesimd_dv stddev_vec(stddev);   

    start = std::chrono::high_resolution_clock::now();

    for(int i = 0;i < sample_size;i += SIMD_SIZE)
    {
        compute_with_umesimd((umesimd_dv &)(x[i]), (umesimd_dv &)(mean[i]),
                             (umesimd_dv &)(stddev[i]), (umesimd_dv &)(result[i]));
    }

    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Finished UME::SIMD computation at " << std::ctime(&end_time)
              << "Elapsed Time : " << elapsed_seconds.count() << "s\n";
              

    //Naive Scalar computation.
    double result_scalar[sample_size];
    
    start = std::chrono::high_resolution_clock::now();
    compute_with_naive_scalar(x, mean, stddev, result_scalar, sample_size);
    end = std::chrono::high_resolution_clock::now();

    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Finished naive Scalar computation at " << std::ctime(&end_time)
              << "Elapsed Time : " << elapsed_seconds.count() << "s\n";


    //Vdt computation.
    double result_vdt[sample_size];

    start = std::chrono::high_resolution_clock::now();
    compute_with_vdt(x, mean, stddev, result_vdt, sample_size);
    end = std::chrono::high_resolution_clock::now();

    elapsed_seconds = end - start;
    end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Finished Vdt computation at " << std::ctime(&end_time)
              << "Elapsed Time : " << elapsed_seconds.count() << "s\n";


    return 0;
}
