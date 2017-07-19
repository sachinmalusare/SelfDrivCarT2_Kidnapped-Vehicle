/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  num_particles = 500;
  std::random_device gen;
  std::normal_distribution<double> dx(x,std[0]);
  std::normal_distribution<double> dy(y,std[1]);
  std::normal_distribution<double> dtheta(theta,std[2]);

  for (int i = 0; i<num_particles; i++)
  {
    Particle tPart = {i, dx(gen), dy(gen), dtheta(gen), 1.0}; // id, x, y, theta, weight
    particles.push_back (tPart);
  }
  weights.resize(num_particles);
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  std::random_device gen;
  std::normal_distribution<double> dx(0,std_pos[0]*delta_t);
  std::normal_distribution<double> dy(0,std_pos[1]*delta_t);
  std::normal_distribution<double> dyaw(0,std_pos[2]*delta_t);

  for (int i = 0; i<num_particles; i++)
  {
    Particle tPart = particles[i]; // get particle
    if (yaw_rate == 0)
    {
      tPart.x = tPart.x + velocity*delta_t*cos(tPart.theta);
      tPart.y = tPart.y + velocity*delta_t*sin(tPart.theta);
      tPart.theta = tPart.theta;
    }
    else
    {
      tPart.x = tPart.x + velocity/yaw_rate*(sin(tPart.theta+yaw_rate*delta_t)-sin(tPart.theta)) + dx(gen);
      tPart.y = tPart.y + velocity/yaw_rate*(cos(tPart.theta)-cos(tPart.theta+yaw_rate*delta_t)) + dy(gen);
      tPart.theta = tPart.theta + yaw_rate*delta_t + dyaw(gen);
    }

    particles[i] = tPart;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
  if (observations.size()>0)
  {
    for (int i = 0; i < observations.size(); i++)
    {
      double xdist = (observations[i].x-predicted[0].x);
      double ydist = (observations[i].y-predicted[0].y);
      int min_dist_id = predicted[0].id;
      double min_dist = sqrt(xdist*xdist+ydist*ydist);
      double dist = 0.0;
      for (int j = 1; j<predicted.size(); j++)
      {
        xdist = (observations[i].x-predicted[j].x);
        ydist = (observations[i].y-predicted[j].y);
        dist = sqrt(xdist*xdist+ydist*ydist);
        if (dist < min_dist)
        {
          min_dist = dist;
          min_dist_id = predicted[j].id;
        }
      }
      observations[i].id = min_dist_id;
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
      /**
       * updateWeights Updates the weights for each particle based on the likelihood of the
       *   observed measurements.
       * @param sensor_range Range [m] of sensor
       * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
       *   standard deviation of bearing [rad]]
       * @param observations Vector of landmark observations
       * @param map Map class containing map landmarks
       */

	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

  for (int i = 0; i<num_particles; i++)
  {
    std::vector<LandmarkObs> transfered_map_landmarks;

    Particle tPart = particles[i]; // get particle
    double px = tPart.x;
    double py = tPart.y;
    double ptheta = tPart.theta;

    for (int j = 0; j<map_landmarks.landmark_list.size(); j++)
    {
      Map::single_landmark_s temp_landmark_input = map_landmarks.landmark_list[j];
      LandmarkObs temp_landmark_output;
      temp_landmark_output.id = temp_landmark_input.id_i;
      double sin_ptheta = sin(ptheta - M_PI/2); // angle change according to README
      double cos_ptheta = cos(ptheta - M_PI/2); // angle change according to README
      double delta_x = temp_landmark_input.x_f - px;
      double delta_y = temp_landmark_input.y_f - py;
      temp_landmark_output.x = -delta_x*sin_ptheta + delta_y*cos_ptheta;
      temp_landmark_output.y = -delta_x*cos_ptheta - delta_y*sin_ptheta;
      transfered_map_landmarks.push_back(temp_landmark_output);
    }

    std::sort(transfered_map_landmarks.begin(), transfered_map_landmarks.end());

//    ParticleFilter::dataAssociation(transfered_map_landmarks, observations);

    long double weight = 1.0;
    for(int k = 0; k < observations.size(); k++)
    {
      double obs_x = observations[k].x;
      double obs_y = observations[k].y;
      double landmark_x = transfered_map_landmarks[k].x;
      double landmark_y = transfered_map_landmarks[k].y;
//      int obs_id = observations[k].id;
//      double landmark_x = transfered_map_landmarks[obs_id].x;
//      double landmark_y = transfered_map_landmarks[obs_id].y;

      long double denomin = 1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]);
      long double x_nom = ((landmark_x-obs_x) * (landmark_x-obs_x)) / (std_landmark[0] * std_landmark[0]);
      long double y_nom = ((landmark_y-obs_y) * (landmark_y-obs_y)) / (std_landmark[1] * std_landmark[1]);
      long double expon = exp(-0.5*(x_nom+y_nom));
      long double meas_likelyhood = denomin * expon;

      weight *= meas_likelyhood;
    }
    particles[i].weight = weight;
    weights[i] = weight;
  } // End particle loop
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<Particle> new_particles;

  double max_weight = *max_element(begin(weights), end(weights));

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> ddouble(0.0, 2.0*max_weight);
  std::uniform_int_distribution<int> dint(0, num_particles);

  // draw with some noise
  std::normal_distribution<double> dx(0,0.003);
  std::normal_distribution<double> dy(0,0.003);
  std::normal_distribution<double> dtheta(0,0.0001);

  int idx = dint(gen);
  double beta = 0.0;

  for (int i = 1; i<num_particles; i++)
  {
    beta = beta + ddouble(gen);
    while (beta > particles[idx].weight)
    {
      beta = beta - particles[idx].weight;
      idx = (idx+1)%num_particles;
    }
    Particle particle;
    particle.x = particles[idx].x+dx(gen);
    particle.y = particles[idx].y+dy(gen);
    particle.theta = particles[idx].theta+dtheta(gen);
    //particle.x = particles[idx].x;
    //particle.y = particles[idx].y;
    //particle.theta = particles[idx].theta;

    particle.weight = particles[idx].weight;
    new_particles.push_back(particle);
  }
  particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
