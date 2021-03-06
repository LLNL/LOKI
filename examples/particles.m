function write_particle_data
   fout = 'test.h5';
   num_particles = 4;
   h5create(fout, '/num_particles', [1]);
   h5write(fout, '/num_particles', num_particles);
   particle_data = zeros(7, num_particles);
   particle_data(1, 1) = 0.0;
   particle_data(2, 1) = -1.0;
   particle_data(3, 1) = 1.0;
   particle_data(4, 1) = 5.0;
   particle_data(5, 1) = 104.0;
   particle_data(6, 1) = 2.0;
   particle_data(7, 1) = 5.0;
   particle_data(1, 2) = 0.0;
   particle_data(2, 2) = -1.0;
   particle_data(3, 2) = 1.0;
   particle_data(4, 2) = 7.0;
   particle_data(5, 2) = 120.0;
   particle_data(6, 2) = 3.0;
   particle_data(7, 2) = 2.0;
   particle_data(1, 3) = 0.0;
   particle_data(2, 3) = -1.0;
   particle_data(3, 3) = 1.0;
   particle_data(4, 3) = 8.0;
   particle_data(5, 3) = 86.0;
   particle_data(6, 3) = 2.0;
   particle_data(7, 3) = 8.0;
   particle_data(1, 4) = 0.0;
   particle_data(2, 4) = -1.0;
   particle_data(3, 4) = 1.0;
   particle_data(4, 4) = 3.0;
   particle_data(5, 4) = 95.0;
   particle_data(6, 4) = 4.0;
   particle_data(7, 4) = 3.0;
   h5create(fout, '/particle_data', [7, 4]);
   h5write(fout, '/particle_data', particle_data);
end
