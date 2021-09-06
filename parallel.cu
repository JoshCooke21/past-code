#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

const double pi = 4.0 * atan(1.0);
__constant__ double c_pi;
__constant__ double c_spacing;
__constant__ int c_num_dipoles;

// Define this to turn on error checking
#define CUDA_ERROR_CHECK

// define wrapper functions giving file and line number of error
#define cuda_safe_call( err ) __cuda_safe_call( err, __FILE__, __LINE__ )
#define cuda_check_error() __cuda_check_error( __FILE__, __LINE__ )

inline void __cuda_safe_call(cudaError err,const char *file,const int line){
  // only produce code if error checking enabled
  #ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err ){ // check for error
    // print out error message
    std::cerr << "cuda_safe_call() failed at " << file << ":" << line
    << " " << cudaGetErrorString( err ) << std::endl;
    // exit code with error code -1
    exit( -1 );
  }
  #endif
  return;
}

class vtimer_t{

private:
   std::chrono::high_resolution_clock::time_point start_time;
   std::chrono::high_resolution_clock::time_point end_time;

public:
   // start the timer
   void start(){
      start_time = std::chrono::high_resolution_clock::now();
   }

   // start the timer
   void stop(){
      end_time = std::chrono::high_resolution_clock::now();
   }

   // get the elapsed time in milliseconds
   double elapsed_time(){

      // work out elapsed time
      return 1.e-9*double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

   }
};



__global__
void kernel_function(double *dev_x_coords, double *dev_y_coords, double *dev_z_coords,
                     double *dev_unit_mu, double *dev_b_x, double *dev_b_y, double *dev_b_z){

  //This array will hold the B field contributions from every dipole onto A per thread
  // to then be summed on threadIdx=0 and entered into the dev_b_(x,y,z) vector,
  // which holds the B field value at every dipole
  extern __shared__ double shared_b[];

  //The outer for loop will be governed by the number of blocks
  //The inner for loop                     the number of threads
  int outer_start  = blockIdx.x;
  int outer_stride = gridDim.x;

  int inner_start  = threadIdx.x;
  int inner_stride = blockDim.x;

  //Pull constants into local memory
  int d_num_dipoles = c_num_dipoles;
  double d_spacing = c_spacing;
  double d_pi = c_pi;

  if (threadIdx.x==0 && blockIdx.x==0) printf("Running with %d blocks and %d threads/block \n", gridDim.x, blockDim.x);

  for (int i=outer_start; i<d_num_dipoles; i+=outer_stride){ // Loop over dipole A

    double coords_i[3]; //Stores the position of dipole A
    coords_i[0] = dev_x_coords[i];
    coords_i[1] = dev_y_coords[i];
    coords_i[2] = dev_z_coords[i];

    //Each thread has its own version of addition that it keeps adding to
    //  so that the sum of addition across all threads in a block gives
    //  the magnetic field at that dipole
    double addition[3]={};

    //Calculate the self-interaction on tid=0 so it only gets added once
    if(threadIdx.x==0){
      double spacing3 = d_spacing*d_spacing*d_spacing;
      for (int n=0; n<3; n++) addition[n] += 2.0 * dev_unit_mu[n] / (3.0 * spacing3 );
    }

    for (int j=inner_start; j<d_num_dipoles; j+=inner_stride){ // Loop over dipole B

      if (i!=j){ // Non-self Interaction

        double disp_vec[3] = {dev_x_coords[j] - coords_i[0], dev_y_coords[j] - coords_i[1], dev_z_coords[j] - coords_i[2]};

        double vec_length = sqrt(disp_vec[0]*disp_vec[0] + disp_vec[1]*disp_vec[1] + disp_vec[2]*disp_vec[2]);
        double vec_length3 = vec_length*vec_length*vec_length;

        double unit_disp_vec[3]={};
        for (int n=0; n<3; n++){ unit_disp_vec[n] = disp_vec[n] / vec_length; }

        double dot_product = unit_disp_vec[0]*dev_unit_mu[0] + unit_disp_vec[1]*dev_unit_mu[1] + unit_disp_vec[2]*dev_unit_mu[2];

        for (int n=0; n<3; n++){
          addition[n] += 0.25 * (1.0/d_pi) * (3.0 * unit_disp_vec[n] * dot_product - dev_unit_mu[n] ) / (vec_length3);
        }

      } // Self-interaction is already accounted for on threadId 0

    } // Close dipole B loop

    //Ensure all threads have finished loop 2
    __syncthreads();

    //Put the contributions from each thread into the shared array
    shared_b[threadIdx.x] = addition[0];
    shared_b[threadIdx.x+blockDim.x] = addition[1];
    shared_b[threadIdx.x+2*blockDim.x] = addition[2];

    if (threadIdx.x==0){

      //Sum up the contributions from all thread calculations
      double dipole_mag_x = 0.0;
      double dipole_mag_y = 0.0;
      double dipole_mag_z = 0.0;
      for (int n=0; n<blockDim.x; n++) dipole_mag_x += shared_b[n];
      for (int n=0; n<blockDim.x; n++) dipole_mag_y += shared_b[blockDim.x+n];
      for (int n=0; n<blockDim.x; n++) dipole_mag_z += shared_b[2*blockDim.x+n];

      //Store the magnetic field at dipole A
      dev_b_x[i] = dipole_mag_x;
      dev_b_y[i] = dipole_mag_y;
      dev_b_z[i] = dipole_mag_z;

    } // Close if threadId==0

  } // Close dipole A loop

  //Ensure all threads have finished all calculations
  __syncthreads();
}


int main(int argc, char* argv[]){

  vtimer_t timer;
  timer.start();

  std::vector<double> unit_mag_moment = {0.0 , 0.0 , 1.0}; // Set magnetic moment direction, currently in z_hat
  std::vector<double> total_mag(3, 0.0); // Initialize the total magnetic field of the ovoid

  // Ellipsoid parameters
  const double r_x=100.0, r_y=100.0, r_z=200.0; // Ellipsoid radii, 10nm and 20nm
  const double spacing=3; // Spacing of points (3 Angstroms)

  //These arrays store the coordinates of the dipoles inside the ellipsoid
  std::vector<double> x_coords(0), y_coords(0), z_coords(0);

  // Looping over z first means that the centre-plane dipoles are grouped together in the coordinate vectors
  for (int k=-ceil(r_z/spacing); k<(ceil(r_z/spacing)+1); k++){ // Looping over dipoles in z
    double z = double(k)*spacing;
    for (int i=-ceil(r_x/spacing); i<(ceil(r_x/spacing)+1); i++){ // Looping over dipoles in x
      double x = double(i)*spacing;
      for (int j=-ceil(r_y/spacing); j<(ceil(r_y/spacing)+1); j++){ // Looping over dipoles in y
        double y = double(j)*spacing;

        if ( x*x/(r_x*r_x) + y*y/(r_y*r_y) + z*z/(r_z*r_z) <= 1){
          x_coords.push_back(x);
          y_coords.push_back(y);
          z_coords.push_back(z);
        }

      } // Close k
    } // Close j
  } // Close i

  //All three coordinate arrays have the same size
  int num_dipoles = x_coords.size();
  std::cout << "Number of dipoles inside: " << num_dipoles << std::endl;
  std::cout << std::endl;

  //These vectors will store the values of the magnetic field at each dipole
  std::vector<double> b_x_alldipoles(num_dipoles, 0.0);
  std::vector<double> b_y_alldipoles(num_dipoles, 0.0);
  std::vector<double> b_z_alldipoles(num_dipoles, 0.0);

  //Allocate arrays on device
  printf("Allocating arrays on device\n");
  int bytes = num_dipoles*sizeof(double); //Same for all three x,y,z arrays
  double *dev_x_coords, *dev_y_coords, *dev_z_coords, *dev_unit_mu, *dev_b_x, *dev_b_y, *dev_b_z;
  cuda_safe_call(cudaMalloc(&dev_x_coords, bytes));
  cuda_safe_call(cudaMalloc(&dev_y_coords, bytes));
  cuda_safe_call(cudaMalloc(&dev_z_coords, bytes));
  cuda_safe_call(cudaMalloc(&dev_unit_mu, 3*sizeof(double)));
  cuda_safe_call(cudaMalloc(&dev_b_x, bytes));
  cuda_safe_call(cudaMalloc(&dev_b_y, bytes));
  cuda_safe_call(cudaMalloc(&dev_b_z, bytes));

  //Copy data from host to device
  printf("Copying data to device\n");
  cuda_safe_call(cudaMemcpy(dev_x_coords, x_coords.data(), bytes, cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_y_coords, y_coords.data(), bytes, cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_z_coords, z_coords.data(), bytes, cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_unit_mu, unit_mag_moment.data(), 3*sizeof(double), cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_b_x, b_x_alldipoles.data(), bytes, cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_b_y, b_y_alldipoles.data(), bytes, cudaMemcpyHostToDevice));
  cuda_safe_call(cudaMemcpy(dev_b_z, b_z_alldipoles.data(), bytes, cudaMemcpyHostToDevice));

  cuda_safe_call(cudaMemcpyToSymbol(c_spacing, &spacing, sizeof(double)));
  cuda_safe_call(cudaMemcpyToSymbol(c_num_dipoles, &num_dipoles, sizeof(int)));
  cuda_safe_call(cudaMemcpyToSymbol(c_pi, &pi, sizeof(double)));

  //Run the parallelized section
  int num_threads_in_block = std::stoi(argv[1]);
  int num_blocks = std::stoi(argv[2]);

  //This sets the number of bytes for the shared array
  int total_shared_bytes = 3*num_threads_in_block*sizeof(double);

  std::cout << "Beginning calculation." << std::endl;

  kernel_function<<<num_blocks, num_threads_in_block, total_shared_bytes>>>(dev_x_coords, dev_y_coords, dev_z_coords, dev_unit_mu,
                                                                            dev_b_x, dev_b_y, dev_b_z);

  //Ensure the CPU and GPU are at the same point (GPU finished, CPU waiting)
  cudaDeviceSynchronize();
  std::cout << "Calculation complete." << std::endl << std::endl;

  //Copy data from device to host
  cuda_safe_call(cudaMemcpy(&b_x_alldipoles[0], dev_b_x, bytes, cudaMemcpyDeviceToHost));
  cuda_safe_call(cudaMemcpy(&b_y_alldipoles[0], dev_b_y, bytes, cudaMemcpyDeviceToHost));
  cuda_safe_call(cudaMemcpy(&b_z_alldipoles[0], dev_b_z, bytes, cudaMemcpyDeviceToHost));

  printf("Copied data back from device \n");

  std::ofstream dfile;
  dfile.open("centre_plane_p.dat");
  dfile << "# x, y, B_z" << std::endl;

  //These are needed to return B to SI units
  const double mu   = 9.274*1E-24;
  const double mu_0 = 4.0*pi*1E-7;

  //Since all the z=0 points are grouped, I use new_x and old_x to check when the x-coordinate changes
  //This allows the 2D plot to plot correctly
  printf("Writing to centre_plane.dat\n");
  double old_x = 0;
  for (int i=0; i<num_dipoles; i++){
    if (z_coords[i]==0){
      double new_x = x_coords[i];

      if (old_x!=new_x) dfile << std::endl; // The newline allows the plotting
      // If on centre plane of ellipsoid, write x and y coords and z-component of dipole_mag
      dfile << x_coords[i] << " " << y_coords[i] << " " << b_z_alldipoles[i]* mu*mu_0*1E30 << std::endl;

      old_x = new_x;
    }
  }

  dfile.close();
  printf("Finished writing to centre_plane.dat\n");

  printf("Totalling the magnetic field values\n");
  for (int n=0; n<num_dipoles; n++){
    total_mag[0]+=b_x_alldipoles[n];
    total_mag[1]+=b_y_alldipoles[n];
    total_mag[2]+=b_z_alldipoles[n];
  }

    // Calculate average magnetisation in all 3 directions (store as vector)
  printf("Calculating the demagnetizing factor\n\n");
  std::vector<double> demagnet(3,0.0);
  std::cout << "Average B, Demagnet. factor" << std::endl;
  
  for (int n=0; n<3; n++){
    double average_mag = total_mag[n] / x_coords.size();
    demagnet[n] = 1.0 - pow(spacing,3) * average_mag;
    std::cout << average_mag* mu*mu_0*1E30 << " " << demagnet[n] << std::endl;;
  }
  std::cout << std::endl; // Prints all vector elements on the same line


  cudaFree(dev_x_coords);
  cudaFree(dev_y_coords);
  cudaFree(dev_z_coords);
  cudaFree(dev_unit_mu);
  cudaFree(dev_b_x);
  cudaFree(dev_b_y);
  cudaFree(dev_b_z);

  timer.stop();
  // printf("%d %d %f \n", num_blocks, num_threads_in_block, timer.elapsed_time());
  std::cout << "Calculation complete. Total elapsed time: " << timer.elapsed_time() << std::endl;

  return 0;

}
