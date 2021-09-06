#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

const double pi = 4.0 * atan(1.0);
const double mu   = 9.274*1E-24;
const double mu_0 = 4.0*pi*1E-7;

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

int main(){

  vtimer_t timer;
  timer.start();

  std::vector<double> unit_mag_moment = {0.0 , 0.0 , 1.0}; // Set magnetic moment direction, currently in z_hat
  std::vector<double> total_mag(3, 0.0); // Set the total magnetic field of the ovoid

  // Ellipsoid parameters
  const double r_x=100.0, r_y=100.0, r_z=200.0; // Ellipsoid radii, 10nm and 20nm
  const double spacing=3; // Spacing of points (in Angstroms)

  std::vector<double> x_coords(0), y_coords(0), z_coords(0);

  // Looping over z first means that the centre-plane dipoles are grouped together in the coordinate vectors
  for (int k=-ceil(r_z/spacing); k<(ceil(r_z/spacing)+1); k++){ // Looping over dipoles in z
    double z = double(k)*spacing;
    for (int i=-ceil(r_x/spacing); i<(ceil(r_x/spacing)+1); i++){ // Looping over dipoles in x
      double x = double(i)*spacing;
      for (int j=-ceil(r_y/spacing); j<(ceil(r_y/spacing)+1); j++){ // Looping over dipoles in y
        double y = double(j)*spacing;

        if ( pow( x/r_x ,2) + pow( y/r_y ,2) + pow( z/r_z ,2) <= 1){
          x_coords.push_back(x);
          y_coords.push_back(y);
          z_coords.push_back(z);
        }

      } // Close k
    } // Close j
  } // Close i

  std::cout << "Number of dipoles inside: " << x_coords.size() << std::endl;
  std::cout << std::endl;
  std::cout << "Beginning calculation." << std::endl;

  //This vector of vectors stores the values of the magnetic field at every dipole
  std::vector< std::vector<double> > b_alldipoles(0);

  for (int i=0; i<x_coords.size(); i++){ // Loop over dipole A

    std::vector<double> dipole_mag(3,0.0);
    std::vector<double> coords_i = {x_coords[i], y_coords[i], z_coords[i]};

    //Add the self-interaction
    double spacing3 = spacing*spacing*spacing;
    for (int n=0; n<3; n++) dipole_mag[n] += 2.0 * unit_mag_moment[n] / (3.0 * spacing3 );

    for (int j=0; j<x_coords.size(); j++){ // Loop over dipole B

      if (i!=j){ // Non-self interaction

        std::vector<double> disp_vec = {x_coords[j] - coords_i[0], y_coords[j] - coords_i[1], z_coords[j] - coords_i[2]};

        double vec_length = sqrt(disp_vec[0]*disp_vec[0] + disp_vec[1]*disp_vec[1] + disp_vec[2]*disp_vec[2]);
        double vec_length3 = vec_length*vec_length*vec_length;

        std::vector<double> unit_disp_vec(3,0.0);
        for (int n=0; n<3; n++){ unit_disp_vec[n] = disp_vec[n] / vec_length; }

        double dot_product = unit_disp_vec[0]*unit_mag_moment[0] + unit_disp_vec[1]*unit_mag_moment[1] + unit_disp_vec[2]*unit_mag_moment[2];

        for (int n=0; n<3; n++){
          dipole_mag[n] += 0.25 * (1.0/pi) * (3.0 * unit_disp_vec[n] * dot_product - unit_mag_moment[n] ) / (vec_length3);
        }
      }

    } // Close dipole B FOR

    //Store the value of the magnetic field
    b_alldipoles.push_back(dipole_mag);

    //Add the magnetic field values to a total
    for (int n=0; n<3; n++) total_mag[n] += dipole_mag[n];
  }

  std::ofstream dfile;
  dfile.open("centre_plane_s.dat");
  dfile << "# x, y, B_z" << std::endl;

  //Since all the z=0 points are grouped, I use new_x and old_x to check when the x-coordinate changes
  //This allows the 2D plot to plot correctly
  double old_x = 0;
  for (int i=0; i<x_coords.size(); i++){
    if (z_coords[i]==0){
      double new_x = x_coords[i];

      if (old_x!=new_x) dfile << std::endl; // The newline allows the plotting
      // If on centre plane of ellipsoid, write x and y coords and z-component of dipole_mag
      dfile << x_coords[i] << " " << y_coords[i] << " " << b_alldipoles[i][2] * mu*mu_0*1E30 << std::endl;

      old_x = new_x;
    }
  }

  dfile.close();

  // Calculate average magnetisation in all 3 directions (store as vector)
  std::vector<double> demagnet(3,0.0);
  std::cout << "Average B, Demagnet. factor" << std::endl;

  for (int n=0; n<3; n++){
    double average_mag = total_mag[n] / x_coords.size();
    demagnet[n] = 1.0 - spacing*spacing*spacing * average_mag;
    std::cout << average_mag * mu*mu_0*1E30 << " " << demagnet[n] << std::endl;;
    // The factor of mu*mu_0*1E30 puts the field in SI
  }

  timer.stop();
  std::cout << "Program complete. Elapsed time: " << timer.elapsed_time() << std::endl;

  return 0;

}
